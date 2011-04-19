classdef BasePCDSystem < models.BaseDynSystem & ISimConstants
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
    % @todo Set typical concentrations at model level (state scaling)
    
    properties       
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
    
    properties(Dependent)
        % Spatial stepwidth
        h; % is set in subclasses
    end
    
    properties(Access=private)
        fh = [];
    end
    
    properties(SetAccess=private)
        % Relative diffusion coefficient (d2/d1)
        D2;
        
        % Relative diffusion coefficient (d3/d1)
        D3;
        
        % Relative diffusion coefficient (d4/d1)
        D4;
        
        % scaled spatial stepwidth
        hs;
    end
    
    methods
        function this = BasePCDSystem(model)
            this = this@models.BaseDynSystem(model);
            
            % Scale diffusion coefficients
            m = this.Model;
            this.D2 = m.d2/m.d1;
            this.D3 = m.d3/m.d1;
            this.D4 = m.d4/m.d1;
            
            this.x0 = @(mu)this.initialX(mu);
            
            % Set output conversion
            %this.C = dscomponents.PointerOutputConv(@(t,mu)this.getC(t,mu), false);

            this.updateSimConstants;
%             this.setParam('Kc1', this.Kc1_real, 1); % *this.ya0
%             this.setParam('Kc2', this.Kc2_real, 1); % *this.xa0^this.n
%             this.setParam('Kd1', this.Kd1_real, 1);
%             this.setParam('Kd2', this.Kd2_real, 1);
%             this.setParam('Kd3', this.Kd3_real, 1);
%             this.setParam('Kd4', this.Kd4_real, 1);
%             this.setParam('Kp1', this.Kp1_real, 1); %/this.xi0
%             this.setParam('Kp2', this.Kp2_real, 1); %/this.yi0
        end
        
        function checkCFL(this)
            m = this.Model;
            if ~isa(m.ODESolver,'solvers.ode.MLImplSolver')
                mi = max([m.d1 m.d2 m.d3 m.d4]);
                if mi*m.dt > .95*this.h^2
                    m.dt = .95*this.h^2/mi;
                    fprintf('Attention. CFL condition violated, setting model dt=%5.3e\n',m.dt);
                end
            end
        end
        
        function set.h(this, value)
            % Set and check CFL
            this.fh = value;
            this.checkCFL;
            this.hs = value / this.Model.L;
            this.updateDims;
        end
        
        function h = get.h(this)
            h = this.fh;
        end
        
        function updateSimConstants(this)
            t = this.Model.tau;
            this.setParam('Kc1', this.Kc1_real * t, 1); % *this.ya0
            this.setParam('Kc2', this.Kc2_real * t, 1); % *this.xa0^this.n
            this.setParam('Kd1', this.Kd1_real * t, 1);
            this.setParam('Kd2', this.Kd2_real * t, 1);
            this.setParam('Kd3', this.Kd3_real * t, 1);
            this.setParam('Kd4', this.Kd4_real * t, 1);
            this.setParam('Kp1', this.Kp1_real * t, 1); %/this.xi0
            this.setParam('Kp2', this.Kp2_real * t, 1); %/this.yi0
            
            this.checkCFL;
        end
        
    end
    
    methods(Abstract, Access=protected)
        updateDims;
        
        x0 = initialX(mu);
        
        C = getC(t,mu);
    end

end



