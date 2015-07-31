classdef DynLinTimoshenkoSystem < models.BaseSecondOrderSystem
% DynLinTimoshenkoSystem: 
%
%
%
% @author Daniel Wirtz @date 2011-09-20
%
% @new{0,5,dw,2011-09-20} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing

    properties(SetAccess=private)
        % The stiffness matrix for t=0, used in the linear model and also
        % in the nonlinear model to compute C.
        K0;
    end
        
    methods
        function this = DynLinTimoshenkoSystem(model)
            this = this@models.BaseSecondOrderSystem(model);
            
            this.MaxTimestep = [];
            
            %% Modellparameter
            % d1 = 0.3 D�mpfungsfaktor vor Massenmatrix (Luftwiderstand)
            this.addParam('d1',1,'Range',[0 2]);
            % d2 = .01 D�mpfungsfaktor vor Steifigkeitsmatrix (Materiald�mpfung)
            this.addParam('d2',.05,'Range',[0 .1]);
            % Add the 3rd parameter, the local gravity factor
            this.addParam('gravity constant',9.81,'Range',[-20 20]);
            
            % Input coordinate one is for constant force term b (previously in Ax+b in corefun)
            % Fake constant gravity
            this.Inputs{1} = @(t)[1; 0; 0; -1];
            this.Inputs{2} = @(t)[1; -1; 0; 0];
            this.Inputs{3} = @(t)[1; 0; -1; 0];
            % More "advanced" gravity
            this.Inputs{4} = @(t)[1; sin(t); 0; 0];
            this.Inputs{5} = @(t)[1; sin((t/10)*2*pi); cos((t/10)*2*pi); 0];
            this.Inputs{6} = @(t)[1; 0; sin((t/10)*2*pi); cos((t/10)*2*pi)];
            
            this.buildElementDependentComponents;
        end
        
        function buildElementDependentComponents(this)
            % Creates all beam element-dependent components like mass matrix and force
            % conversion matrix.
            m = this.Model;
            data = m.data;
            
            %% Assemblieren der Massenmatrix & Steifigkeitsmatrix f�r t=0
            % Mass matrix (M is actually the full M (including dirichlet
            % nodes), but is projected at the end of this method to the
            % size of actual DoFs)
            e = m.Elements;
            i = []; j = [];
            M = []; K = [];
            for k = 1:length(e)
                M_lok = e{k}.getLocalMassMatrix;
                K_lok = e{k}.getLocalStiffnessMatrix;
                index_glob = e{k}.getGlobalIndices;

                [li,lj] = meshgrid(index_glob);
                i = [i; li(:)]; %#ok<*AGROW>
                j = [j; lj(:)];
                M = [M; M_lok(:)];
                K = [K; K_lok(:)];
            end
            nk = data.num_knots;
            
            %% Mass matrix
            M = sparse(i,j,M,7*nk,7*nk);
            % Project M to effectively needed entries
            M = M(m.free,m.free);
            this.M = dscomponents.ConstMassMatrix(M);
            % Create full K0 matrix
            this.K0 = sparse(i,j,K,7*nk,7*nk);
            K = this.K0(m.free,m.free);
            
            %% Stiffness matrix
            if m.NonlinearModel
                this.f = models.beam.DLTNonlinearCoreFun(this);
            else
                this.A = dscomponents.LinearCoreFun(-K);
            end
            
            %% Fill inner AffLinCoreFun with matrices
            % Dämpfungsmodell 1: M a_t + (d1*M + d2*K) v_t + K u_t = f
            % d1 = 0.3 Dämpfungsfaktor vor Massenmatrix (Luftwiderstand)
            % d2 = .01 Dämpfungsfaktor vor Steifigkeitsmatrix (Materiald�mpfung)
            d = dscomponents.AffLinCoreFun(this);
            d.addMatrix('mu(1)',M);
            d.addMatrix('mu(2)',-K);
            this.D = d;
            
            dim = length(m.free);
            this.NumStateDofs = dim;
            this.NumDerivativeDofs = dim;
            this.updateDimensions;
            
            %% Initial values
            this.x0 = dscomponents.ConstInitialValue(zeros(dim,1));
                        
            %% Model input (gravity + const term)
            % Affine-parametric matrix models const. 
            % B(t,\mu) = \mu_3*[0 F] + 1*[f_const 0 0 0], u(t) = [1 g(t)], g\in\R^3
            B = dscomponents.AffLinInputConv;
            
            %% 1*[f_const 0 0 0]
            % Neumann forces (computed in Model.preprocess_data)
            f_const = m.f_neum;
            f_const = f_const(m.free);
            % Dirichlet forces only for linear model
            if ~m.NonlinearModel
                % Reconstruct fake u vector which has entries only at
                % dirichlet points. Used for computation of constant forces due
                % to dirichlet values.
                u = zeros(7*nk,1);
                u(m.dir_u) = m.u_dir;
                f_dir = -this.K0*u;
                f_const = f_const + f_dir(m.free);
            end
            Fconst = sparse(dim,4);
            Fconst(:,1) = f_const;
            B.addMatrix('1',Fconst);
            
            %% \mu_3*[0 F]
            i = []; j = [];
            F = []; 
            for k = 1:length(e)
                F_lok = e{k}.getLocalForceMatrix';
                index_glob = e{k}.getGlobalIndices;

                [li,lj] = meshgrid(index_glob,1:3);
                i = [i; li(:)]; %#ok<*AGROW>
                j = [j; lj(:)];
                F = [F; F_lok(:)];
            end
            
            F = sparse(i,j,F,7*nk,3);
            F = F(m.free,:);
            %F = [sparse(size(F,1),3); F];
            % Add zero column to front
            F = [sparse(size(F,1),1) F];
            B.addMatrix('mu(3)',F);
            B.CoeffClass = 'InputConvCoeffs';
            
            this.B = B;
            
            %% Output setup (dim = 3*space + 3*velo + 1*temp)
            %xdim = (dim - (dim/7))/2;
            nf = length(m.free);
            C = speye(nf);
            % Remove velocity entries
            %C(repmat(logical([0 0 0 1 1 1]'),nf/6,1),:) = [];
            C = [C 0*C];
            this.C = dscomponents.LinearOutputConv(C);
            % Export settings
            je = m.JKerMorExport;
            je.DoFFields = 6;
            
            f = struct;
            f(1).Type = 'Displacement3D';
            f(1).Name = [];
            f(2).Type = 'RealValue';
            f(2).Name = 'X-Velocity';
            f(3).Type = 'RealValue';
            f(3).Name = 'Y-Velocity';
            f(4).Type = 'RealValue';
            f(4).Name = 'Z-Velocity';
            je.LogicalFields = f;
            
            this.updateSparsityPattern;
        end
    end
    
    methods(Access=protected)
        function val = getDerivativeDirichletValues(this, t)
            val = [];
        end
    end
end