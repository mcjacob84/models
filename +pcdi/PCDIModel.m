classdef PCDIModel < models.BaseFullModel
    % Base class inhibitor PCD models
    %
    % @author Daniel Wirtz @date 2013-10-23
    %
    % @new{0,8,dw,2013-11-22} Moved SteadyStates into a function that
    % computes the steady states depending on the current parameter
    % (exponent mu(4))
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
            
            if nargin < 2
                inhibitors = true;
                if nargin < 1
                    dim = 2;
                end
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
        
        function ss = getSteadyStates(this, mu)
            % Computes the steady state configurations depending on the
            % current system (given by mu)
            %
            % Note: These are SCALED quantities. The cell is considered
            % dead as of ya > 0.01.
            %
            % Parameters:
            % mu: The current parameter @type colvec<double>
            %
            % Return values:
            % ss: The steady states. First to third row: life state,
            % unstable state, death state
            
            % Hack: bar seems to be a script file that accidentally gets
            % called even though its initialized as variable.. so fake
            % assign stuff to make the parser think its a variable.
            bar = [];
            syms xa ya xi yi iap bar yb xb;
            
            % Reaction system
            F = [this.K2*xi*ya-this.K5*xa-this.K11*xa*bar+this.Km11*xb;...
                this.K1*yi*xa^mu(4)-this.K6*ya-this.K3*ya*iap+this.Km3*yb;
                -this.K2*xi*ya-this.K9*xi+this.Km9;...
                -this.K1*yi*xa^mu(4)-this.K10*yi+this.Km10;...
                -this.K3*ya*iap-this.K8*iap+this.Km8-this.K4*ya*iap+this.Km3*yb;...
                -this.K11*xa*bar+this.Km11*xb-this.K12*bar+this.Km12;...
                this.K3*ya*iap-this.Km3*yb-this.K7*yb;...
                this.K11*xa*bar-this.Km11*xb-this.K13*xb];
            sol = solve(F);
            % Find solution with real and positive concentrations
            pos = find((sol.xa >= 0) .* (sol.ya >= 0) .* (sol.xi >= 0) .* ...
                (sol.yi >= 0) .* (sol.iap >= 0) .* (sol.bar >= 0) .* ...
                (sol.yb >= 0) .* (sol.xb >= 0));
            if length(pos) == 3
                ss = [[sol.xa(pos(1));sol.xa(pos(2));sol.xa(pos(3))]...
                    [sol.ya(pos(1));sol.ya(pos(2));sol.ya(pos(3))]...
                    [sol.xi(pos(1));sol.xi(pos(2));sol.xi(pos(3))]...
                    [sol.yi(pos(1));sol.yi(pos(2));sol.yi(pos(3))]...
                    [sol.iap(pos(1));sol.iap(pos(2));sol.iap(pos(3))]...
                    [sol.bar(pos(1));sol.bar(pos(2));sol.bar(pos(3))]...
                    [sol.yb(pos(1));sol.yb(pos(2));sol.yb(pos(3))]...
                    [sol.xb(pos(1));sol.xb(pos(2));sol.xb(pos(3))]];
            else
                % Life State independent of the exponent n
                ss = [0 0 this.Km9/this.K9 this.Km10/this.K10...
                    this.Km8/this.K8 this.Km12/this.K12 0 0];
            end
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

