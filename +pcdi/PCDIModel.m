classdef PCDIModel < models.BaseFullModel
    % Base class inhibitor PCD models
    %
    % @author Daniel Wirtz @date 2013-10-23
    %
    % @new{0,8,dw,2013-10-23} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    
    properties(Constant)
        %% System Coefficient values
        % Typical diffusion rate for Caspase-8
        d1 = 1.8e-11; %[m^2/s]
        
        % Typical diffusion rate for Caspase-3
        d2 = 1.86e-11; %[m^2/s]
        
        % Typical diffusion rate for Pro-Caspase-8
        d3 = 1.89e-11; %[m^2/s]
        
        % Typical diffusion rate for Pro-Caspase-3
        d4 = 2.27e-11; %[m^2/s]
        
        % Typical diffusion rate for IAP
        d5 = 1.8720e-11; %[m^2/s]
        
        % Typical diffusion rate for BAR
        d6 = 1.9548e-11; %[m^2/s]
        
        % Typical diffusion rate for YAI
        d7 = 1.4850e-11; %[m^2/s]
        
        % Typical diffusion rate for XAB
        d8 = 1.4850e-11; %[m^2/s]
        
        % Procaspase-8 to Caspase-8 reaction rate
        K1 = 0.08;
        
        % Procaspase-3 to Caspase-3 reaction rate
        K2 = 0.08; 
        
        % IAP-Caspase 3 (de)reaction rate
        K3 = 0.054; 
        
        % IAP-Caspase 3 one-way (de)reaction rate
        K4 = 0.36; 
        
        % Caspase-8 degradation rate
        K5 = .0005; 
        
        % Caspase-3 degradation rate
        K6 = .0005; 
        
        % YAI degradation rate
        K7 = 3.6000e-4;
        
        % IAP degradation rate
        K8 = 1.8000e-04;
        
        % Pro-Caspase-8 degradation rate
        K9 = .0005;
        
        % Caspase-3 degradation rate
        K10 = .0005;
        
        % BAR - Procaspase-8 (de)reaction rate
        K11 = 0.54; 
        
        % BAR degradation rate
        K12 = 1.8000e-05; 
        
        % XAP degradation rate
        K13 = 1.8000e-4;
        
        % YAI to IAP production rate
        Km3 = 0.0360;
        
        % IAP production rate
        Km8 = 1e-04; 
        
        % Procaspase-8 production rate
        Km9 = 1e-04; 
        
        % Procaspase-3 production rate
        Km10 = 1e-04;
        
        % XAB degradation rate
        Km11 = 0.0360;
        
        % BAR production rate
        Km12 = 1.8000e-06;
        
        %% System Rescaling settings
        % Typical concentration
        tc = 1e-7; %[Mol]
        
        % Typical cell length (from 1D)
        L = 1e-5; %[m]
        
        % Steady state configurations
        %
        % Note: These are SCALED quantities
        %
        % First row: life state
        % Second row: unstable state
        % Third row: death state (as of ya > 0.01 its considered death)
        SteadyStates = [[0; 2.16e-4; 0.0097]...
            [0; 0.0012; 0.1109]...
            [.1984; 0.1659; 0.0107]...
            [.1984; 0.1918; 0.0781]...
            [5.5556; 2.4177; 0.0473]...
            [1; 0.8372; 0.1027]...
            [0; 0.0015; 0.0026]...
            [0; 0.009; 0.0498]]
        % Determination of steady states for the given parameter values
        % dependent on the exponent mu(4)
        %
        % syms xa ya xi yi iap bar yb xb;
        % Reaction system
        % F = [K2*xi*ya-K5*xa-K11*xa*bar+Km11*xb;...
        %     K1*yi*xa^mu(4)-K6*ya-K3*ya*iap+Km3*yb;
        %     -K2*xi*ya-K9*xi+Km9;...
        %     -K1*yi*xa^mu(4)-K10*yi+Km10;...
        %     -K3*ya*iap-K8*iap+Km8-K4*ya*iap+Km3*yb;...
        %     -K11*xa*bar+Km11*xb-K12*bar+Km12;...
        %     K3*ya*iap-Km3*yb-K7*yb;...
        %     K11*xa*bar-Km11*xb-K13*xb];
        % Sol = solve(F);
        
        % Find solution with real and positive concentrations
        % Index = find((Sol.xa >= 0) .* (Sol.ya >= 0) .* (Sol.xi >= 0) .* ...
        % (Sol.yi >= 0) .* (Sol.iap >= 0) .* (Sol.bar >= 0) .* ...
        % (Sol.yb >= 0) .* (Sol.xb >= 0))
        % if length(Index) == 3
        %    SteadyStates =
        %    [[Sol.xa(Index(1));Sol.xa(Index(2));Sol.xa(Index(3))]...
        %     [Sol.ya(Index(1));Sol.ya(Index(2));Sol.ya(Index(3))]...
        %     [Sol.xi(Index(1));Sol.xi(Index(2));Sol.xi(Index(3))]...
        %     [Sol.yi(Index(1));Sol.yi(Index(2));Sol.yi(Index(3))]...
        %     [Sol.iap(Index(1));Sol.iap(Index(2));Sol.iap(Index(3))]...
        %     [Sol.bar(Index(1));Sol.bar(Index(2));Sol.bar(Index(3))]...
        %     [Sol.yb(Index(1));Sol.yb(Index(2));Sol.yb(Index(3))]...
        %     [Sol.xb(Index(1));Sol.xb(Index(2));Sol.xb(Index(3))]...
        % else
        % Life State independent of the exponent n
        %     SteadyStates = [0 0 Km9/K9 Km10/K10 Km8/K8 Km12/K12 0 0];
        % end                                
    end
    
    properties(SetAccess=private)
        WithInhibitors;
    end
    
    properties(Dependent)
        Dimension;
    end
    
    methods
        function this = PCDIModel(dim, inhibitors)
            % Creates a new instance of the PCDModel
            %
            % Parameters:
            % dim: The dimension to use @default 1
            
            if nargin == 0
                dim = 2;
            end
            this.WithInhibitors = inhibitors;
            
            this.T = 3600; %[s]
            this.dt = 5; %[s]
            % time scaling
            this.tau = this.L^2/this.d1;
            
            s = sampling.GridSampler;
            this.Sampler = s;
            
            %             this.ODESolver = solvers.MLWrapper(@ode23);
            %             this.ODESolver = solvers.ExplEuler(this.dt);
            o = solvers.MLode15i;
            o.AbsTol = 1e-6;
            o.RelTol = 1e-5;
            o.MaxStep = [];
            this.ODESolver = o;
            %             this.ODESolver = solvers.SemiImplicitEuler(this);
            
            switch dim
                case 3
                    s = models.pcdi.PCDISystem3D(this);
                    this.Name = 'Programmed Cell Death 3D';
                    this.Data.useFileTrajectoryData;
                otherwise
                    s = models.pcdi.PCDISystem2D(this);
                    this.Name = 'Programmed Cell Death 2D';
                    this.Data.useFileTrajectoryData;
            end
            if this.WithInhibitors
                this.Name = [this.Name ' (Inhibitor System)'];
            end            
            
            s.MaxTimestep = this.dt;
            this.System = s;
            
            % Space reduction setup
            sr = spacereduction.PODGreedy;
            sr.Eps = 1e-5;
            this.SpaceReducer = sr;
            
            this.Approx = [];
        end
        
        function plot(this, varargin)
            % Overrides standard method and forwards to the system's plot
            % function. (they are 1D and 2D)
            this.System.plot(this, varargin{:});
        end
        
        function plotState(this, varargin)
            % Overrides standard method and forwards to the system's plot
            % function. (they are 1D and 2D)
            this.System.plotState(this, varargin{:});
        end
        
        function value = get.Dimension(this)
            value = this.System.Dims;
        end
        
        function set.WithInhibitors(this, value)
            if ~islogical(value) || ~isscalar(value)
                error('WithInhibitors must be a logical scalar');
            end
            this.WithInhibitors = value;
        end
    end
end

