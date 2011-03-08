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
    
    properties       
        % The parent PCD model
        model;
        
        % Spatial stepwidth
        h = .1;
        
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
        % Typical cell length (from 1D)
        L = 1e-5; %[m]
        
        % Typical diffusion rate for Caspase-8
        d1 = 1.8e-11; %[m^2/s]
        % Typical diffusion rate for Caspase-3
        d2 = 1.86e-11; %[m^2/s]
        % Typical diffusion rate for Pro-Caspase-8
        d3 = 1.89e-11; %[m^2/s]
        % Typical diffusion rate for Pro-Caspase-3
        d4 = 2.27e-11; %[m^2/s]
    end
    
    properties(SetAccess=private)        
        % Relative diffusion coefficient (d2/d1)
        D2;
        
        % Relative diffusion coefficient (d3/d1)
        D3;
        
        % Relative diffusion coefficient (d4/d1)
        D4;
    end
    
    methods
        function this = BasePCDSystem(model)
            this.model = model;
            this.x0 = @(mu)this.initialX(mu);
            
            % Set output conversion
            %this.C = dscomponents.PointerOutputConv(@(t,mu)this.getC(t,mu), false);
            
            this.setParam('Kc1', this.Kc1_real, 1); % *this.ya0
            this.setParam('Kc2', this.Kc2_real, 1); % *this.xa0^this.n
            this.setParam('Kd1', this.Kd1_real, 1);
            this.setParam('Kd2', this.Kd2_real, 1);
            this.setParam('Kd3', this.Kd3_real, 1);
            this.setParam('Kd4', this.Kd4_real, 1);
            this.setParam('Kp1', this.Kp1_real, 1); %/this.xi0
            this.setParam('Kp2', this.Kp2_real, 1); %/this.yi0
        end
        
        function updateSimConstants(this)
            % Initializes constants for a simulation.
            %this.lam1 = this.xi0/this.xa0;
            %this.lam2 = this.yi0/this.ya0;
            t = this.L^2/this.d1;
            % Set scaling of the model from here
            this.model.tau = t;
            this.D2 = this.d2/this.d1;
            this.D3 = this.d3/this.d1;
            this.D4 = this.d4/this.d1;
            this.updateDims;
            
            this.setParam('Kc1', this.Kc1_real * t, 1); % *this.ya0
            this.setParam('Kc2', this.Kc2_real * t, 1); % *this.xa0^this.n
            this.setParam('Kd1', this.Kd1_real * t, 1);
            this.setParam('Kd2', this.Kd2_real * t, 1);
            this.setParam('Kd3', this.Kd3_real * t, 1);
            this.setParam('Kd4', this.Kd4_real * t, 1);
            this.setParam('Kp1', this.Kp1_real * t, 1); %/this.xi0
            this.setParam('Kp2', this.Kp2_real * t, 1); %/this.yi0
        end
        
    end
    
    methods(Abstract, Access=protected)
        updateDims;
        x0 = initialX(mu);
        C = getC(t,mu);
    end
end



