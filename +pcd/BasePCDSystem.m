classdef BasePCDSystem < models.BaseDynSystem
    %PCDSYSTEM The 2D dynamical system of the Programmed Cell Death Model
    %by Markus Daub.
    %
    %   For details on properties and components see the Paper
    %   "Mathematical Modeling of Programmed Cell Death", Markus Daub,
    %   Milestone presentation.
    %
    % See also: models.pcd.CoreFun
    %
    % @author Daniel Wirtz @date 2010-03-15
    %
    % @new{0,6,dw,2011-11-27} New property BasePCDSystem.SteadyStates that includes the
    % linearized system's steady states for life and death (From "M. Daub: Death wins against
    % life in a spatially extended apoptosis model"). This property is used in the plot methods
    % for the different models.
    %
    % @new{0,5,dw,2011-11-02} 
    % - New property ReacCoeff. The coefficients have been removed (commented out) of the
    % system's Params and are inserted manually into the mu parameters in the CoreFunctions in
    % order not to have to change much code. This way, they can be easily reintroduced as true
    % system parameters.
    % - Removed checkCFL method and placed it into the set.h method (directly sets the
    % MaxTimestep property according to CFL conditions)
    %
    % @change{0,5,dw,2011-10-10} Generalized and restructured the models.pcd models and
    % systems as far as possible. Removed the ISimConstants interfaces and
    % added protected, custom models.pcd.BasePCDSystem.newSysDimension methods to propagate changes
    % in geometry (models.pcd.BasePCDSystem.Omega or
    % models.pcd.BasePCDSystem.h) to the systems and core functions.
    %
    % @new{0,3,dw,2011-04-21} Integrated this class to the property default value changed
    % supervision system @ref propclasses. This class now inherits from KerMorObject and has an
    % extended constructor registering any user-relevant properties using
    % KerMorObject.registerProps.
    
    properties(SetObservable, Dependent)
        % Spatial stepwidth (in unscaled size units!)
        %
        % @propclass{critical} Determines the spatial resolution of the model.
        %
        % See also: L
        h; % is set in subclasses
        
        % The spatial width/area/region (in unscaled size units!)
        %
        % @propclass{data} The area for the model.
        %
        % See also: L
        Omega;
    end
    
    properties(Access=private)
        fh = [];
        fOmega = [];
    end
    
    properties(SetAccess=private)
        % Relative diffusion coefficients ([d2/d1, d3/d1, d4/d1])
        Diff;
        
        % scaled spatial stepwidth
        hs;
        
        % The system's dimensions
        Dims;
        
        % The reaction coefficients
        ReacCoeff;
    end
    
    methods
        function this = BasePCDSystem(model)
            this = this@models.BaseDynSystem(model);
            
            % Scale diffusion coefficients
            m = this.Model;
            this.Diff = [m.d2/m.d1 m.d3/m.d1 m.d4/m.d1];
            
            this.registerProps('h','Omega');
            
            this.ReacCoeff = [m.Kc1_real m.Kc2_real m.Kd1_real ...
                              m.Kd2_real m.Kd3_real m.Kd4_real ...
                              m.Kp1_real m.Kp2_real]' * m.tau;
                          
            % Add parameters
            % Activation area
            p = this.addParam('area', [.01, 1], 10);
            p.Spacing = 'lin';
            
            % Activation rate
            p = this.addParam('rate', [0.0005, 1], 15);
            p.Spacing = 'log';
            
            % Activation time
            p = this.addParam('atime', [0, 1], 5);
            p.Spacing = 'lin';
        end
        
        function h = get.h(this)
            h = this.fh;
        end
        
        function set.h(this, value)
            if any(value >= this.fOmega(:,2))
                error('Cannot choose a step size h=%e value larger or equal to the geometry [%s].',...
                    value,Utils.implode(this.fOmega(:),', ','%g'));
            end
            this.fh = value;
            this.hs = value / this.Model.L;
            this.updateDims;
            
            m = this.Model;
            this.MaxTimestep = [];
            if ~isa(m.ODESolver,'solvers.AImplSolver') && ~isa(m.ODESolver,'solvers.SemiImplicitEuler')
                maxdt = .95*(value^2/max([m.d1 m.d2 m.d3 m.d4]));
                maxdtsc = .95*(this.hs^2/max([1 this.Diff]));
                
                if (maxdtsc - maxdt/this.Model.tau) /  maxdtsc > 1e-15
                    error('Inconsistent scaling. Please check.');
                end
                % Set max timestep value (scaled!)
                this.MaxTimestep = maxdtsc;
                
                if maxdt < m.dt
                    warning('models:pcd:CFL','CFL condition violated with h=%e and current dt=%e.\nSetting System.MaxTimestep=%e (scaled, effective value %e)\n', value, m.dt, maxdtsc, maxdt);
                end
            end
        end
        
        function v = get.Omega(this)
            v = this.fOmega;
        end
        
        function set.Omega(this, value)
            if size(value,2) ~= 2
                error('Omega needs to be a dim x 2 matrix with the spatial extend for each dimension in a row.');
            end
            this.fOmega = value;
            
            % Update the dimensions
            this.updateDims;
        end
        
        function setConfig(this, mu, inputidx)
            setConfig@models.BaseDynSystem(this, mu, inputidx);
            m = this.Model;
            this.ReacCoeff = [m.Kc1_real m.Kc2_real m.Kd1_real ...
                              m.Kd2_real m.Kd3_real m.Kd4_real ...
                              m.Kp1_real m.Kp2_real]' * m.tau;
        end
    end
    
    methods(Access=private)
        
        function updateDims(this)
            if ~isempty(this.fh) && ~isempty(this.fOmega)
                nd = size(this.fOmega,1);
                this.Dims = zeros(1,nd);
                for d=1:nd
                    this.Dims(d) = length(this.fOmega(d,1):this.h:this.fOmega(d,2));
                end
                m = prod(this.Dims);
                
                % Set state scaling
                ss = zeros(4*m,1);
                ss(1:m) = this.Model.xa0;
                ss(m+1:2*m) = this.Model.ya0;
                ss(2*m+1:3*m) = this.Model.xi0;
                ss(3*m+1:end) = this.Model.yi0;
                this.StateScaling = ss;
                
                % Set new initial values and output in 1D-3D systems
                this.newSysDimension;
                
                % Call template method for custom actions in the core
                % functions (diffusion matrix comp)
                this.f.newSysDimension;
            end
        end
    end
    
    methods(Abstract, Access=protected)
        % Custom updates for new system dimension
        newSysDimension(this);
    end

end



