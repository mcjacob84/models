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
        % Spatial stepwidth
        h = .1;
        
        % Exponent in ya,yi term (necessary casp-3 for casp-8 activation)
        n = 2;
        
        %% Coefficient values
        
        % Procaspase-8 to Caspase-8 reaction rate
        % Empiric value from [1] in daub's milestone
        Kc1_real = 1e5;
        
        % Procaspase-3 to Caspase-3 reaction rate
        % Empiric value from [1] in daub's milestone
        Kc2_real = 5.8e12;
        
        % Caspase-8 degradation rate
        % Empiric value from [1] in daub's milestone
        Kd1_real = .0001;
        
        % Caspase-3 degradation rate
        % Empiric value from [1] in daub's milestone
        Kd2_real = .0001;
        
        % Procaspase-8 production rate
        % Empiric value from [1] in daub's milestone
        Kp1_real = 1.4e-11;
        
        % Procaspase-3 production rate
        % Empiric value from [1] in daub's milestone
        Kp2_real = 2.3e-12;
        
        %% System Rescaling settings
        
        % Typical Caspase-8 concentration
        xa0 = 1e-7; %[M]
        % Typical Caspase-3 concentration
        ya0 = 1e-7; %[M]
        % Typical Procaspase-8 concentration
        xi0 = 1e-7; %[M]
        % Typical Procaspase-3 concentration
        yi0 = 1e-7; %[M]
        % Typical Length
        L = 1e-5; %[m]
        % Typical diffusion rate for (Pro-)Caspase-8
        d1 = 1e-11; %[m^2/s]
        % Typical diffusion rate for (Pro-)Caspase-3
        d2 = 1e-11; %[m^2/s]
    end
    
    properties(SetAccess=private)
        % Procaspase-8 to Caspase-8 typical concentration quotient
        lam1;
        % Procaspase-3 to Caspase-3 typical concentration quotient
        lam2;
        % Typical timestep (dependent)
        tau; %[s]
        % Relative diffusion coefficient (d1/d2)
        D;
    end
    
    methods
        function this = BasePCDSystem
            this.x0 = @(mu)this.initialX(mu);
            
            % Set output conversion
            this.C = dscomponents.PointerOutputConv(@(t,mu)this.getC(t,mu), false);
        end
        
        function updateSimConstants(this)
            % Initializes constants for a simulation.
            this.lam1 = this.xi0/this.xa0;
            this.lam2 = this.yi0/this.ya0;
            this.tau = this.L^2/this.d1;
            this.D = this.d1/this.d2;
            this.updateDimSimConstants;
            
            this.setParam('Kc1', this.Kc1_real * this.ya0*this.tau, 1);
            this.setParam('Kc2', this.Kc2_real * this.xa0^this.n*this.tau, 1);
            this.setParam('Kd1', this.Kd1_real * this.tau, 1);
            this.setParam('Kd2', this.Kd2_real * this.tau, 1);
            this.setParam('Kp1', this.Kp1_real * this.tau/this.xi0, 1);
            this.setParam('Kp2', this.Kp2_real * this.tau/this.yi0, 1);
        end
        
    end
    
    methods(Abstract, Access=protected)
        updateDimSimConstants;
        x0 = initialX(mu);
        C = getC(t,mu);
    end
end



