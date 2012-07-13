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
    %
    % @todo Set typical concentrations at model level (state scaling)
    
    properties(Constant)
        % Exponent in ya,yi term (necessary casp-3 for casp-8 activation)
        n = 2;
        
        %% Coefficient values        
        % Procaspase-8 to Caspase-8 reaction rate
        % Empiric value from [1] in daub's milestone
        Kc1_real = 0.08;
        
        % Procaspase-3 to Caspase-3 reaction rate
        % Empiric value from [1] in daub's milestone
        Kc2_real = 0.08;
        
        % Caspase-8 degradation rate
        % Empiric value from [1] in daub's milestone
        Kd1_real = .0005;
        
        % Caspase-3 degradation rate
        % Empiric value from [1] in daub's milestone
        Kd2_real = .0005;
        
        % Pro-Caspase-8 degradation rate
        % Empiric value from [1] in daub's milestone
        Kd3_real = .0005;
        
        % Caspase-3 degradation rate
        % Empiric value from [1] in daub's milestone
        Kd4_real = .0005;
        
        % Procaspase-8 production rate
        % Empiric value from [1] in daub's milestone
        Kp1_real = 0.0001;
        
        % Procaspase-3 production rate
        % Empiric value from [1] in daub's milestone
        Kp2_real = 0.0001;
        
        %% System Rescaling settings
        % Typical Caspase-8 concentration
        xa0 = 1e-7; %[M]
        
        % Typical Caspase-3 concentration
        ya0 = 1e-7; %[M]
        
        % Typical Procaspase-8 concentration
        xi0 = 1e-7; %[M]
        
        % Typical Procaspase-3 concentration
        yi0 = 1e-7; %[M]
        
        % Steady state configurations
        % First row: life state
        % Second row: unstable state
        % Third row: death state
        SteadyStates = [[0; 9.8153e-4; 0.1930]*models.pcd.BasePCDSystem.xa0...
                        [0; 3.0824e-5; 0.1713]*models.pcd.BasePCDSystem.ya0...
                        [.2; 0.1990; 0.0070]*models.pcd.BasePCDSystem.xi0...
                        [.2; 0.2; 0.0287]*models.pcd.BasePCDSystem.yi0];
    end
    
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
            
            % Uncomment if reaction parameters become real parameters again
            %this.setReactionParams;
            this.ReacCoeff = [this.Kc1_real this.Kc2_real this.Kd1_real ...
                              this.Kd2_real this.Kd3_real this.Kd4_real ...
                              this.Kp1_real this.Kp2_real]' * this.Model.tau;
        end
        
        function h = get.h(this)
            h = this.fh;
        end
        
        function set.h(this, value)
            if any(value >= this.fOmega(:,2))
                error('Cannot choose a step size h=%e value larger or equal to the geometry [%s].',value,num2str(this.fOmega));
            end
            this.fh = value;
            this.hs = value / this.Model.L;
            this.updateDims;
            
            m = this.Model;
            this.MaxTimestep = [];
            if ~isa(m.ODESolver,'solvers.ode.AImplSolver') && ~isa(m.ODESolver,'solvers.ode.SemiImplicitEuler')
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
            
            %this.setReactionParams;
            this.ReacCoeff = [this.Kc1_real this.Kc2_real this.Kd1_real ...
                              this.Kd2_real this.Kd3_real this.Kd4_real ...
                              this.Kp1_real this.Kp2_real]' * this.Model.tau;
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
                ss(1:m) = this.xa0;
                ss(m+1:2*m) = this.ya0;
                ss(2*m+1:3*m) = this.xi0;
                ss(3*m+1:end) = this.yi0;
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



