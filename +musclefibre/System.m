classdef System < models.motorunit.MotorunitBaseSystem
% System: The global dynamical system used within the Model
%
% Contains the Dynamics, Linear Diffusion and InputConv components as well as a
% ConstInitialValue.
%
% @author Daniel Wirtz @date 2012-11-22
%
% @new{0,8,dw,2015-09-15} Imported into +models package.
%
% @new{0,7,dw,2012-11-22} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    properties(SetAccess=private)
        % The number of sarcomer cells in this fibre model
        % 
        % @type integer @default 100
        % See also: musclefibres.Model
        N;
        
        % The number of the sarcomere that gets the input of the motoneuron
        %
        % @type integer
        MotoSarcoLinkIndex;
        
        ds = 0; % No spindle as of now
    end
    
    properties(Dependent)
        % Length of the muscle fibre [cm]
        %
        % Computed from N*dx
        %
        % See also: N dx
        FibreLength;
    end
    
    properties(Constant) 
        % Diffusion coefficient
        % conductivity tensor (single entry in 1D) in [mS/cm]
        sigma = 3.828;
        
        % surface to volume ratio of a muscle fibre in [cm^-1]
        % (assume all fibres have same length and diameter)
        A_m = 500;
        
        % "length" of one sarcomere unit [cm]
        %
        % The default value is oriented on Thomas Heidlauf's work.
        %
        % @type double @default 0.0052
        dx = 6/1152; % [cm]
    end
    
    properties(Access=private)
        spindle;
    end
    
    methods
        function this = System(model, options)
            % Call superclass constructor
            this = this@models.motorunit.MotorunitBaseSystem(model, options);
            
            % Set local variables
            this.N = options.N;
            
            % neuromuscular junction makes more sense (symmetry...)
            %this.MotoSarcoLinkIndex = round(options.N/2);
            this.MotoSarcoLinkIndex = 1;
            
            this.NumStateDofs = this.dsa*options.N+this.dm+this.ds;
            this.updateDimensions;
            
            % Linear diffusion part
            this.assembleA;
            
            %% Set system components
            % Core nonlinearity
            this.f = models.musclefibre.Dynamics(this, options);
            
            % Linear input B for motoneuron
            this.assembleB;
            
            % Affine-Linear output C
            this.assembleC;
            
            % Initial values
            this.assembleX0(options);
            
            this.updateSparsityPattern;
        end
        
        function value = get.FibreLength(this)
            value = this.N*this.dx;
        end
    end
    
    methods(Access=protected)
        
        function A = assembleA(this)
            % Computes the linear term A which represents the diffusion of
            % the membrane voltage along the sarcomeres. Diffusion is
            % discretized via finite differences. Since membrane voltage is
            % one of 58 states, A is very sparse. Diffusion depends on fibre type (fast - slow). 
            
            n0=this.dm+this.ds;  % size of upper left zeros block
            sar=this.dsa; % dimension of single sarcomer part
            size=this.dm+this.ds+sar*this.N;

            i=[n0+1+sar:sar:size,n0+1:sar:size,n0+1:sar:size-sar];    % row index
            j=[n0+1:sar:size-sar,n0+1:sar:size,n0+1+sar:sar:size];   % column index
            s=[ones(1,this.N-1),-1,-2*ones(1,this.N-2),-1,ones(1,this.N-1)];   % values
            s = (this.sigma/this.A_m) * s / this.dx^2;
            As = sparse(i,j,s,size,size);
            
            
            % Pre-multiply constant part
            A = dscomponents.AffLinCoreFun(this);
            A.TimeDependent = false;
            coeff_str = sprintf('1./(mu(1,:)*%g + (1-mu(1,:))*%g)',...
                this.sarco.C_m_fast,this.sarco.C_m_slow);
            A.addMatrix(coeff_str,As);
            this.A = A;
        end
        
        function B = assembleB(this)
            % input conversion matrix, depends on fibre type. Has only one
            % entry in second row.
            B = dscomponents.AffLinInputConv;
            % Base noise input mapping
            B.addMatrix('1./(pi*(exp(log(100)*mu(1,:))*3.55e-05 + 77.5e-4).^2)',...
                sparse(2,1,1,this.dm+this.ds+this.dsa*this.N,2));
            % Independent noise input mapping with µ_2 as mean current factor
            B.addMatrix('mu(2,:)./(pi*(exp(log(100)*mu(1,:))*3.55e-05 + 77.5e-4).^2)',...
                sparse(2,2,1,this.dm+this.ds+this.dsa*this.N,2));
            this.B = B;
        end
        
        function assembleC(this)
%             % input conversion matrix, depends on fibre type. Has only one
%             % entry in second row.
%             C = dscomponents.AffLinOutputConv;
%             % Extract V_s
%             C.addMatrix('1',sparse(1,7,1,2,this.dm+this.dsa));
%             % Add scaling for force output A_s
%             str = '1';
%             if this.SingleTwitchOutputForceScaling
%                 str = ['polyval([' sprintf('%g ',this.ForceOutputScalingPolyCoeff) '],mu(1))'];
%             end
%             C.addMatrix(str, sparse(2,59,1,2,this.dm+this.dsa));
        end
        
        function assembleX0(this,options)
            sx0 = [];
            if options.Spindle
                sx0 = this.spindle.InitialValues;
            end
            if this.DynamicInitialConditions
                mc = ?models.motorunit.Shorten;
                s = load(fullfile(fileparts(which(mc.Name)),...
                    sprintf('x0coeff%d.mat',this.SarcoVersion)));
                x0 = dscomponents.AffineInitialValue;
                m = this.dm+this.ds+this.N*this.dsa;
                % The first this.dm coeffs are for the motoneuron
                for k=1:this.dm
                    x0.addMatrix(sprintf('polyval([%s],mu(1))',...
                        sprintf('%g ',s.coeff(k,:))),full(sparse(k,1,1,m,1)));
                end
                % The next this.dsa coeffs are for the sarcomere.
                % Simply repeat them for all sarcomeres!
                for k=1:this.dsa
                    %mat = full(sparse(k,1,1,m,1));
                    pos = this.dm + (k:this.dsa:this.N*this.dsa);
                    mat = sparse(pos,1,1,m,1);
                    x0.addMatrix(sprintf('polyval([%s],mu(1))',...
                        sprintf('%g ',s.coeff(this.dm+k,:))),mat);
                end
                this.x0 = x0;
            else
                this.x0 = dscomponents.ConstInitialValue(...
                    [this.moto.InitialValues;
                    sx0;
                    repmat(this.sarco.InitialValues,this.N,1)]);
            end
        end
        
    end
    
end