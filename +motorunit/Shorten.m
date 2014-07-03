classdef Shorten < models.BaseFullModel
    % Shorten: Model for a muscle motor unit composed of motoneuron
    % and a sarcomere.
    %
    % The global time unit for this model is milliseconds [ms].
    % This model has been copied from the MuscleFibreModel in order to
    % obtain a current snapshot without changing the multi-sarcomere
    % adoptions done for the whole fibre model.
    %
    % This model is intended for research regarding the
    % kernel-approximation of the motor unit w.r.t. the activation and
    % fibre type relation to the force development.
    %
    % @author Daniel Wirtz @date 2014-01-16
    %
    % @new{0,7,dw,2014-01-16} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    
    properties(Dependent)
        UseNoise;
    end
    
    properties(SetAccess=private)
        % Flag that determines if this system is to be run so that only
        % ever one peak/signal will be issued.
        %
        % @type logical @default false
        SinglePeakMode = false;
        
        % Set this flag to true (in constructor) if you want
        % parameter-dependent initial conditions that have been obtained
        % using long-time simulations of the model for different fibre
        % types and then using the end state as "stable" initial condition.
        %
        % The actual values have been computed using a polynomial fit of
        % degree seven on each dimension for the parameter range [0,1]
        %
        % @type logical @default true
        DynamicInitialConditions = true;
    end
    
    methods
        function this = Shorten(dynamic_ic, singlepeakmode)
            % Creates a new motor unit model
            %
            % Parameters:
            % singlepeakmode: A flag that determines if this model should
            % only ever produce one force peak @type logical @default false
            
            if nargin < 2
                singlepeakmode = false;
                if nargin < 1
                    dynamic_ic = true;
                end
            end
            this.SinglePeakMode = singlepeakmode;
            this.DynamicInitialConditions = dynamic_ic;
            
            if singlepeakmode
                this.T = 150; % [ms]
            else
                this.T = 500; % [ms]
            end
            % DO NOT INCREASE! Peaks from Motoneuron are not correctly
            % resolved (at least visually) if larger timesteps are used.
            this.dt = .1; % [ms]
            
            this.SaveTag = 'motorunit';
            this.Data = data.ModelData(this);
            this.Data.useFileTrajectoryData;
            
            this.Name = 'Motor unit model';
            this.System = models.motorunit.SHSystem(this);
            this.TrainingInputs = 1;
            
            s = solvers.MLWrapper(@ode15s);
            if singlepeakmode
                sys = this.System;
                s.odeopts.OutputFcn = @sys.singlePeakModeOutputFcn;
            end
            this.ODESolver = s;
            
            % Set parameter domain already
            s = sampling.RandomSampler;
            s.Domain = models.motoneuron.ParamDomain;
            this.Sampler = s;
            
            this.DefaultMu = [.1; 3];
            this.DefaultInput = 1;
        end
        
        function pm = plotMotoSacroLinkFactorCurve(this)
            x = 0:.1:80;
            pm = PlotManager;
            pm.LeaveOpen = true;
            h = pm.nextPlot('moto_sarco_link_factor','Factor for motoneuro to sarcomere link','Moto V_s','Factor');
            f = this.System.f;
            fx = f.MSLink_MaxFactor*ones(1,length(x));
            dynfac = x < f.MSLink_MaxFactorSignal;
            fx(dynfac) = f.getLinkFactor(x(dynfac));
            plot(h,x,fx);
            pm.done;
        end
        
        function plotOutputForceScaling(this, x)
            if nargin < 2
                x = 0:.01:1;
            end
            plot(x,polyval(this.System.ForceOutputScalingPolyCoeff,x));
            title('Force scaling curve for different fibre types');
            xlabel('Fibre type parameter');
            ylabel('Peak force for single excitation');
        end
        
        function pm = plotState(this, t, x, pm)
            if nargin < 4
                pm = PlotManager(false,3,1);
                pm.LeaveOpen = true;
            end
            h = pm.nextPlot('moto','Motoneuron V_s','time','value');
            plot(h,t,x(2,:));
            h = pm.nextPlot('sarco','Linked sarcomere: V_s','time','V_s');
            plot(h,t,x(this.System.dm+1,:));
            h = pm.nextPlot('sarco',sprintf('Linked sarcomere: A_2\nMu=[%s]',num2str(this.System.mu')),'time','A_2');
            plot(h,t,x(this.System.dm+53,:));
            
            if nargin < 4
                pm.done;
            end
        end
        
        function pm = plot(this, t, y, pm)
            if nargin < 4
                pm = PlotManager(false,2,1);
                pm.LeaveOpen = true;
            end
            h = pm.nextPlot('sarco','Linked sarcomere: V_s','time','V_s');
            plot(h,t,y(1,:));
            h = pm.nextPlot('sarco',sprintf('Linked sarcomere: A_2\nMu=[%s]',num2str(this.System.mu')),'time','A_2');
            plot(h,t,y(2,:));
            
            if nargin < 4
                pm.done;
            end
        end
    end
    
    methods
        function value = get.UseNoise(this)
            value = ~this.System.noiseGen.DisableNoise;
        end
        
        function set.UseNoise(this, value)
            this.System.noiseGen.DisableNoise = ~value;
        end
    end
    
    methods(Static, Access=protected)
        function this = loadobj(this)
            if ~isa(this, 'models.motorunit.Shorten')
                sobj = this;
                this = models.motorunit.Shorten;
                this = loadobj@models.BaseFullModel(this, sobj);
            else
                this = loadobj@models.BaseFullModel(this);
            end
        end
    end 
end