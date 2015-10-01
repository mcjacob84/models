classdef MotorunitBaseSystem < models.BaseFirstOrderSystem
% MotorunitBaseSystem: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2015-09-15
%
% @new{0,7,dw,2015-09-15} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/
% - \c Documentation http://www.morepas.org/software/kermor/
% - \c License @ref licensing
    
    properties(Constant)    
        % Coefficients of the polynomial used for force output scaling.
        %
        % The scaling is applied to the resulting forces, so that the force
        % generated on a single peak is the same for every fibre type.
        %
        % Old version with original initial values and 1000ms simulation
        % time:
        % ForceOutputScalingPolyCoeff = 1e5*[0.142559788512870  -1.052260908967481   3.446953882000253  -6.597894914676725   8.180848326826332  -6.875124492949854   3.967796548608296  -1.549518903381187   0.388781656140931  -0.054972337733183   0.002625376146421   0.000246831617001   0.000019186970434];
        %
        % @type rowvec<double>
        %
        % See also: models.motorunit.experiments.SarcoScaling
        ForceOutputScalingPolyCoeff = {
            1e5*[   0.133229113241298  -0.834127708845393   2.314092280509735  -3.754880326472914   3.945406185349044  -2.792486768513515   1.335923148831657  -0.421435990707753   0.082794510979418  -0.008574151547690  -0.000101330255143  0.000228461898458   0.000022092057889];
            1e10*[  -0.005583080638834   0.051746276520123  -0.220624569577431   0.573660315487840  -1.016764880760808   1.300670698247225  -1.240304735397929   0.897946884962748  -0.498044127370431   0.212068292632907  -0.069047146216459  0.017020761843826  -0.003125097502991   0.000417182230585  -0.000039143308385   0.000002464732214  -0.000000098699610   0.000000002932298   0.000000000218357]};
    end

    properties(SetAccess=private)
        % The 
        SarcoVersion;
        
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
        %
        % See also: models.motorunit.experiment.InitialConditions
        DynamicInitialConditions = true;
        
        % Flag to determine if the output force (e.g. calcium
        % concentration) should be re-scaled so that the forces using a
        % single twitch are the same maximal forces (=1) for each fibre
        % type.
        %
        % @type logical @default true
        %
        % See also: models.motorunit.experiments.SarcoScaling 
        SingleTwitchOutputForceScaling = true;
    end

    properties(SetAccess=private)
        sarco;
        moto;
        
        % Dimension of single sarcomer cell part
        dsa;
        % Dimension of single motoneuron part
        dm;
    
        % The NoiseGenerator to obtain the inputs u(t)        
        noiseGen;
    end
    
    properties(Access=protected, Transient)
        peakstart = false;
    end
    
    properties(Access=private)
        % The upper limit polynomial for maximum mean current dependent on
        % the fibre type.
        %
        % This polynomial has been computed using the
        % models.motoneuron.Model (which is also included here, but has
        % been established as single model for speed and exemplatory
        % purposes), where the fibre type and mean activation current along
        % the 60Hz-contour have been used to fit a polynomial that yields
        % the maximum mean input current for any fibre type.
        %
        % See also: models.motoneuron.Model.FibreTypeDepMaxMeanCurrent
        % models.motoneuron.experiments.ParamDomainDetection
        upperlimit_poly;
    end
    
    methods
        function this = MotorunitBaseSystem(model, options)
            % Call superclass constructor
            this = this@models.BaseFirstOrderSystem(model);
            
            this.SinglePeakMode = options.SPM;
            this.DynamicInitialConditions = options.DynamicIC;
            this.SingleTwitchOutputForceScaling = options.OutputScaling;
            this.SarcoVersion = options.SarcoVersion;
            
            this.addParam('fibre_type', 0, 'Range', [0 1], 'Desired', 15);
            this.addParam('mean_current_factor', 3, 'Range', [0 9], 'Desired', 15);
            
            this.MaxTimestep = model.dt;
            
            % Load mean current limiting polynomial
            s = load(models.motoneuron.Model.FILE_UPPERLIMITPOLY);
            this.upperlimit_poly = s.upperlimit_poly;
            
            % Setup noise input
            ng = models.motoneuron.NoiseGenerator;
            ng.DisableNoise = ~options.Noise;
            this.noiseGen = ng;
            this.Inputs{1} = @ng.getInput;
            
            if options.SarcoVersion == 1
                this.sarco = models.motorunit.SarcomereOriginal;
            else
                this.sarco = models.motorunit.SarcomereNeumann;
            end
            this.dsa = this.sarco.Dims;
            
            this.moto = models.motoneuron.Motoneuron;
            
            this.dm = this.moto.Dims;
        end
        
        function setConfig(this, mu, inputidx)
            % Sets the configuration for the upcoming simulation of the
            % model.
            setConfig@models.BaseFirstOrderSystem(this, mu, inputidx);
             % Precompute motoneuron/sarcomere constants
            this.sarco.setType(mu(1));
            this.moto.setType(mu(1));
            this.noiseGen.setFibreType(mu(1));
            this.MaxTimestep = this.Model.dt*100;
            % Helper method to reset the "peak counter" of the dynamics.
            if this.SinglePeakMode
                this.peakstart = false;
            end
        end
        
        function prepareSimulation(this, mu, inputidx)
            % Limits the mean input current to a maximum value depending on
            % the fibre type, so that the frequency of 60Hz is not
            % exceeded.
            %
            % See also: models.motoneuron.ParamDomainDetection upperlimit_poly
            
            % Limit mean current depending on fibre type
            mu(2) = min(polyval(this.upperlimit_poly,mu(1)),mu(2));
            prepareSimulation@models.BaseFirstOrderSystem(this, mu, inputidx);
        end
        
        function y = computeOutput(this, x, mu)
            if nargin < 3
                mu = this.mu;
            end
            % Subtract the initial value
            x0 = this.getX0(mu);
            x = bsxfun(@minus,x,x0);
            % Get the output
            y = computeOutput@models.BaseFirstOrderSystem(this, x, mu);
            % Set those outputs to zero 
            %y(2,y(2,:)<0) = 0;
        end
        
        function varargout = plot(this, varargin)
            % plots some interesting states of the model
            %
            % See also: musclefibres.System
            [varargout{1:nargout}] = this.System.plot(varargin{:});
        end
    end
    
    methods(Abstract, Access=protected)
        assembleB(this);
        assembleC(this);
        assembleX0(this);
    end
    
end