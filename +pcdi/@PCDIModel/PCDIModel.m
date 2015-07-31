classdef PCDIModel < models.BaseFullModel
    % Base class inhibitor PCD models.
    %
    % The package models.pcdi contains a copy of the models.pcd model
    % files, extended by an optional "inhibitor" extended governing
    % equation. The constructor can be called passing a flag (default=true)
    % if the extension should be used. As this mainly affects the system's
    % nonlinearity, a separate class file has been created to cater for
    % that case. Different behaviour in model and system classes is
    % realized using the flag WithInhibitors.
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
    % - \c Homepage http://www.morepas.org/software/index.html
    % - \c Documentation http://www.morepas.org/software/kermor/index.html
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
        % Flag that detemines if this model uses inhibitors.
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
            % dim: The dimension to use @type integer @default 2
            % inhibitors: Flag to indicate if the inhibitor extensions
            % should be used in this model. @type logical @default true
            
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
    
    methods(Access=protected)
        % Function to determine the steady states, depending on the current
        % parameter mu (i.e. the exponent n which is one parameter)
        %
        % Function defined in separate file in class folder.
        ss = getSteadyStates(this, n);
    end
    
    methods(Static)
        function res = test_PCDIModel_Simulation
            m = models.pcdi.PCDIModel;
            [t,y] = m.simulate;
            m.plot(t,y);
            res = 1;
        end
    end
end

