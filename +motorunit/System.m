classdef System < models.BaseFirstOrderSystem
% SHSystem: The global dynamical system used within the Shorten motorunit
% model
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
    
    properties(Transient)
        % debug output
        d = [];
    end
    
    properties(Constant)
        % Dimension of motoneuron part
        dm = 6; 
        
        % Membrane capacitance. Different values for different fibre types, 
        % due to different action potential propagation speeds
        % C_m is computed parameter-dependent.
        % These constants are for both slow and fast muscle models and are also used in the
        % first entry of the sarcomere constants computed in
        % models.muscle.FibreDynamics.initSarcoConst @type double
        C_m_slow = 0.58;
        C_m_fast = 1;
        
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
        ForceOutputScalingPolyCoeff = 1e5*[0.072141692014344  -0.456658053996476   1.274105896595146  -2.041231359447655   2.058191740934228  -1.353159489096666   0.584711952400098  -0.164458365881107   0.029358733675881  -0.003174883073694   0.000165817182406    0.000048348565452   0.000012193246968];
    end
    
    properties(SetAccess=private)
        % The NoiseGenerator to obtain the inputs u(t)        
        noiseGen;
        
        % Dimension of single sarcomer cell part
        dsa = 56; 
    end
    
    properties(SetAccess=private)
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
    
    properties(Access=private, Transient)
        peakstart = false;
    end
    
    methods
        function this = System(model)
            % Call superclass constructor
            this = this@models.BaseFirstOrderSystem(model);
            
            this.addParam('fibre_type', 0, 'Range', [0 1]);
            this.addParam('mean_current_factor', 3, 'Range', [0 9]);
            
            % The neumann version has two extra equations
            if model.Version == 1
                this.dsa = 56;
            else
                this.dsa = 58;
            end
            
            this.NumStateDofs = this.dm+this.dsa;
            this.updateDimensions;
            
            %% Set system components
            % Core nonlinearity
            this.f = models.motorunit.Dynamics(this);
            
            % Load mean current limiting polynomial
            s = load(models.motoneuron.Model.FILE_UPPERLIMITPOLY);
            this.upperlimit_poly = s.upperlimit_poly;
            
            % Setup noise input
            ng = models.motoneuron.NoiseGenerator;
            this.noiseGen = ng;
            this.Inputs{1} = @ng.getInput;
            
            % Linear input B for motoneuron
            this.B = this.assembleB;
            
            % Affine-Linear output C
            this.C = this.assembleC;
            
            % Constant initial values
            % See also: models.motorunit.experiment.InitialConditions
            if model.DynamicInitialConditions
                this.x0 = this.assembleX0;
            else
                this.x0 = dscomponents.ConstInitialValue(...
                    [this.f.moto.InitialValues;
                    this.f.sarco.InitialValues]);
            end
            this.updateSparsityPattern;
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
        
        function setConfig(this, mu, inputidx)
            % Sets the configuration for the upcoming simulation of the
            % model.
            setConfig@models.BaseFirstOrderSystem(this, mu, inputidx);
            this.noiseGen.setFibreType(mu(1));
            this.MaxTimestep = this.Model.dt*100;
            % Helper method to reset the "peak counter" of the dynamics.
            if this.Model.SinglePeakMode
                this.peakstart = false;
            end
        end
        
        function status = singlePeakModeOutputFcn(this, ~, y, flag)
            if ~strcmp(flag,'done')
                if y(this.dm+1) > -20 && ~this.peakstart
                    this.peakstart = true;
                elseif y(this.dm+1) < -20 && this.peakstart
                    this.peakstart = false;
                    this.f.hadpeak = true;
                end
            end
            status = 0;
        end
    end
    
    methods(Access=private)
        
        function x0 = assembleX0(this)
            % Loads the polynomial coefficients for each dimension and
            % creates an affine-initial value that produces the suitable
            % initial conditions determined by long-time simulations of
            % different fibre-types with no active input.
            %
            % See also: models.motorunit.experiment.InitialConditions
            mc = metaclass(this);
            s = load(fullfile(fileparts(which(mc.Name)),...
                sprintf('x0coeff%d.mat',this.Model.Version)));
            x0 = dscomponents.AffineInitialValue;
            m = size(s.coeff,1);
            for k=1:m
                x0.addMatrix(sprintf('polyval([%s],mu(1))',...
                    sprintf('%g ',s.coeff(k,:))),full(sparse(k,1,1,m,1)));
            end
        end
        
        function B = assembleB(this)
            % input conversion matrix, depends on fibre type. Has only one
            % entry in second row.
            %
            % The divisor in both coefficient functions is the old para.CS
            % value!!
            B = dscomponents.AffLinInputConv;
            % Base noise input mapping
            B.addMatrix('1./(pi*(exp(log(100)*mu(1,:))*3.55e-05 + 77.5e-4).^2)',...
                sparse(2,1,1,this.dm+this.dsa,2));
            % Independent noise input mapping with µ_2 as mean current factor
            B.addMatrix('mu(2,:)./(pi*(exp(log(100)*mu(1,:))*3.55e-05 + 77.5e-4).^2)',...
                sparse(2,2,1,this.dm+this.dsa,2));
        end
        
        function C = assembleC(this)
            % input conversion matrix, depends on fibre type. Has only one
            % entry in second row.
            C = dscomponents.AffLinOutputConv;
            % Extract V_s
            C.addMatrix('1',sparse(1,7,1,2,this.dm+this.dsa));
            % Add scaling for force output A_s
            str = '1';
            if this.Model.SingleTwitchOutputForceScaling
                str = ['polyval([' sprintf('%g ',this.ForceOutputScalingPolyCoeff) '],mu(1))'];
            end
            C.addMatrix(str, sparse(2,59,1,2,this.dm+this.dsa));
        end                
    end
    
end
