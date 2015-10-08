classdef System < models.motorunit.MotorunitBaseSystem
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
    
    methods
        function this = System(model, options)
            % Call superclass constructor
            this = this@models.motorunit.MotorunitBaseSystem(model, options);
            
            this.NumStateDofs = this.dm+this.dsa;
            this.updateDimensions;
            % First (and only) sarcomere index is Vm
            this.MotoSarcoLinkIndex = this.dm +1;
            
            %% Set system components
            % Core nonlinearity
            this.f = models.motorunit.Dynamics(this);
            
            % Linear input B for motoneuron
            this.assembleB;
            
            % Affine-Linear output C
            this.assembleC(options);
            
            % Constant initial values
            this.assembleX0;
            
            this.updateSparsityPattern;
        end
    end
    
    methods(Access=protected)
        
        function assembleX0(this)
            % Loads the polynomial coefficients for each dimension and
            % creates an affine-initial value that produces the suitable
            % initial conditions determined by long-time simulations of
            % different fibre-types with no active input.
            %
            % See also: models.motorunit.experiment.InitialConditions
            
            if this.DynamicInitialConditions
                mc = metaclass(this);
                s = load(fullfile(fileparts(which(mc.Name)),...
                    sprintf('x0coeff%d.mat',this.SarcoVersion)));
                x0 = dscomponents.AffineInitialValue;
                m = size(s.coeff,1);
                for k=1:m
                    x0.addMatrix(sprintf('polyval([%s],mu(1))',...
                        sprintf('%.14g ',s.coeff(k,:))),full(sparse(k,1,1,m,1)));
                end
                this.x0 = x0;
            else
                this.x0 = dscomponents.ConstInitialValue(...
                    [this.moto.InitialValues;
                    this.sarco.InitialValues]);
            end
        end
        
        function assembleB(this)
            % input conversion matrix, depends on fibre type. Has only one
            % entry in second row.
            %
            % The divisor in both coefficient functions is the old para.CS
            % value!!
            B = dscomponents.AffLinInputConv;
            % Base noise input mapping
            B.addMatrix('1./(pi*(exp(log(100)*mu(1,:))*3.55e-05 + 77.5e-4).^2)',...
                sparse(2,1,1,this.NumStateDofs,2));
            % Independent noise input mapping with µ_2 as mean current factor
            B.addMatrix('mu(2,:)./(pi*(exp(log(100)*mu(1,:))*3.55e-05 + 77.5e-4).^2)',...
                sparse(2,2,1,this.NumStateDofs,2));
            this.B = B;
        end
        
        function assembleC(this, options)
            % input conversion matrix, depends on fibre type. Has only one
            % entry in second row.
            C = dscomponents.AffLinOutputConv;
            % Extract V_s
            C.addMatrix('1',sparse(1,7,1,2,this.NumStateDofs));
            % Add scaling for force output A_s
            str = '1';
            coeff = this.ForceOutputScalingPolyCoeff{options.SarcoVersion};
            if this.SingleTwitchOutputForceScaling
                str = ['polyval([' sprintf('%.14g ',coeff) '],mu(1))'];
            end
            C.addMatrix(str, sparse(2,59,1,2,this.NumStateDofs));
            this.C = C;
        end                
    end
    
end
