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
        % Spatial stepwidth
        h = []; % is set in subclasses
        
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
        
        function [bool,mi] = checkCFL(this)
            m = this.Model;
            if ~isa(m.ODESolver,'solvers.MLImplSolver')
                mi = max([m.d1 m.d2 m.d3 m.d4])*m.dt;
                bool = mi < .9*this.h^2;
            else
                bool = true;
                mi = 0;
            end
        end
        
        function set.h(this, value)
            oldh = this.h;
            % Set and check CFL
            this.h = value;
            [b,m] = this.checkCFL;%#ok
            if b
                this.hs = value / this.Model.L;%#ok
                this.updateDims;%#ok
            else
                error('CFL condition violated. Any h²=%5.3e must be greater than %5.3e, consider adjusting the time-step before trying again.',value^2,m);
                this.h = oldh;
            end
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
            
            [b,m] = this.checkCFL;
            if ~b
                fprintf('ATTENTION: CFL condition violated with current dt/h settings. Setting dt = %f\n',.9*m);
                this.Model.dt = .9*m;
            end
        end
        
    end
    
    methods(Abstract, Access=protected)
        updateDims;
        
        x0 = initialX(mu);
        
        C = getC(t,mu);
    end

end



