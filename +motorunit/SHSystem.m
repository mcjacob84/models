classdef SHSystem < models.BaseDynSystem
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
        
        % Dimension of single sarcomer cell part
        dsa = 56; 
        
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
        ForceOutputScalingPolyCoeff = 1e5*[0.072141692014344  -0.456658053996476   1.274105896595146  -2.041231359447655   2.058191740934228  -1.353159489096666   0.584711952400098  -0.164458365881107   0.029358733675881  -0.003174883073694   0.000165817182406    0.000048348565452   0.000012193246968];
    end
    
    properties(SetAccess=private)
        % The NoiseGenerator to obtain the inputs u(t)        
        noiseGen;
    end
    
    properties(Access=private, Transient)
        peakstart = false;
    end
    
    methods
        function this = SHSystem(model)
            % Call superclass constructor
            this = this@models.BaseDynSystem(model);
            
            this.addParam('fibre_type', [0 1], 10);
            this.addParam('mean_current_factor', [.5 10], 10);
            %this.addParam('moto_param', [0 1], 10);
            
            %% Set system components
            % Core nonlinearity
            this.f = models.motorunit.SHDynamics(this);
            
            % Setup noise input
            ng = models.motoneuron.NoiseGenerator;
            this.noiseGen = ng;
            this.Inputs{1} = @ng.getInput;
            
            % Linear input B for motoneuron
            this.B = this.assembleB;
            
            % Affine-Linear output C
            this.C = this.assembleC;
            
            % Constant initial values
            if model.DynamicInitialConditions
                this.x0 = this.assembleX0;
            else
                this.x0 = dscomponents.ConstInitialValue(this.getConstInitialStates);
            end
        end
        
        function setConfig(this, mu, inputidx)
            % Helper method to reset the "peak counter" of the dynamics.
            setConfig@models.BaseDynSystem(this, mu, inputidx);
            this.noiseGen.setMu(mu);
            this.MaxTimestep = this.Model.dt*100;
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
            mc = metaclass(this);
            s = load(fullfile(fileparts(which(mc.Name)),'x0coeff.mat'));
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
            C.addMatrix(['polyval([' sprintf('%g ',this.ForceOutputScalingPolyCoeff) '],mu(1))'],...
                sparse(2,59,1,2,this.dm+this.dsa));
        end
                
        function x0 = getConstInitialStates(~)
            % initial states
            x0_sarc = zeros(58,1);
            x0_sarc(1:14) = [-79.974; -80.2; 5.9; 150.9; 5.9; 12.7; 133;...
                133; 0.009466; 0.9952; 0.0358; 0.4981; 0.581; 0.009466];
            x0_sarc(15:37) = [0.9952; 0.0358; 0.4981; 0.581; 0; 0; 0; 0; 0; 1;...
                0; 0; 0; 0; 0.1; 1500; 0.1; 1500; 25; 615; 615; 811; 811];
            x0_sarc(38:45) = [16900; 16900; 0.4; 0.4; 7200; 7200; 799.6; 799.6];
            x0_sarc(46:54) = [1000; 1000; 3; 0.8; 1.2; 3; 0.3; 0.23; 0.23];
            x0_sarc(55:58) = [0.23; 0.23; 0.0; 0.05];
            x0 = [zeros(6,1); x0_sarc];
            
%             x0 = 1e4*[0.000002105216220
%                0.000003621385002
%                0.000002909625708
%                0.000099590346266
%                0.000003791919729
%                0.000000015081509
%               -0.008220103158330
%               -0.008213826304369
%                0.000616422489154
%                0.015090303609523
%                0.000588455048963
%                0.001269826128002
%                0.013186084039038
%                0.013301169961923
%                0.000000688872253
%                0.000099542136254
%                0.000002706918371
%                0.000059557827555
%                0.000059087570477
%                0.000000695103658
%                0.000099638277984
%                0.000002728603130
%                0.000059289421401
%                0.000059082262288
%                0.000000000176541
%                0.000000000559251
%                0.000000000664353
%                0.000000000350759
%                0.000000000069447
%                0.000088270373229
%                0.000011185010563
%                0.000000531482664
%                0.000000011224301
%                0.000000000088892
%                0.000081729161448
%                0.149193283547150
%                0.000036457584391
%                0.149272225466980
%                0.002645213006267
%                0.061500000000000
%                0.061500000000000
%                0.081100000000000
%                0.081100000000000
%                1.686302951356341
%                1.686709766870293
%                0.000206227516784
%                0.000138893833514
%                0.724232562918369
%                0.724234577079512
%                0.075561209564849
%                0.075626529086977
%                0.095765941567035
%                0.095765438026699
%                0.000395940593721
%                0.000102616321003
%                0.000162655158102
%                0.000391688669685
%                0.000040793527447
%                0.000032519348632
%                0.000023060913139
%                0.000032393965109
%                0.000013164858422
%                                0
%                0.000005000000000];
        end
        
    end
    
end
