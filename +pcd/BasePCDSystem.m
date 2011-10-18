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
    % @author Daniel Wirtz @date 15.03.2010
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
    end
    
    properties(SetObservable, Dependent)
        % Spatial stepwidth
        %
        % @propclass{critical} Determines the spatial resolution of the model.
        h; % is set in subclasses
        
        % The spatial width/area/region
        %
        % @propclass{data} The area for the model.
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
    end
    
    methods
        function this = BasePCDSystem(model)
            this = this@models.BaseDynSystem(model);
            
            % Scale diffusion coefficients
            m = this.Model;
            this.Diff = [m.d2/m.d1 m.d3/m.d1 m.d4/m.d1];
            
            this.registerProps('h','Omega');
            
            this.setReactionParams;
        end
        
        function [res, maxdt] = checkCFL(this, hvalue)
            res = true;
            m = this.Model;
            maxdt = m.dt;
            if ~isa(m.ODESolver,'solvers.ode.MLode15i')
                mi = max([m.d1 m.d2 m.d3 m.d4]);
                if mi*m.dt > .95*hvalue^2
                    maxdt = .95*hvalue^2/mi;             
                    res = false;
                end
            end
        end
        
        function h = get.h(this)
            h = this.fh;
        end
        
        function set.h(this, value)
            % Set and check CFL
            [r,mdt] = this.checkCFL(value);           
            if r
                this.fh = value;
                this.hs = value / this.Model.L;
                this.updateDims;
            else
                error('CFL condition violated, not changing h. Suggested model dt for h=%5.3e is %5.3e\n', value, mdt);
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
            if ~this.checkCFL(this.fh);
                error('CFL condition violated. Check first.');
            end
            
            setConfig@models.BaseDynSystem(this, mu, inputidx);
            
            this.setReactionParams;
        end
    end
    
    methods(Access=private)
        function setReactionParams(this)
            t = this.Model.tau;
            this.setParam('Kc1', this.Kc1_real * t, 1); % *this.ya0
            this.setParam('Kc2', this.Kc2_real * t, 1); % *this.xa0^this.n
            this.setParam('Kd1', this.Kd1_real * t, 1);
            this.setParam('Kd2', this.Kd2_real * t, 1);
            this.setParam('Kd3', this.Kd3_real * t, 1);
            this.setParam('Kd4', this.Kd4_real * t, 1);
            this.setParam('Kp1', this.Kp1_real * t, 1); %/this.xi0
            this.setParam('Kp2', this.Kp2_real * t, 1); %/this.yi0
        end
        
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
        newSysDimension;
    end

end



