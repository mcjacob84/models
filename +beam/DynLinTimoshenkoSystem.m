classdef DynLinTimoshenkoSystem < models.BaseDynSystem
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
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    properties(SetAccess=private)
        % The small mass matrix of half dimension; corresponds to the lower
        % right half block of this.M
        M_small;
        
        % The stiffness matrix for t=0, used in the linear model and also
        % in the nonlinear model to compute C.
        K0;
    end
        
    methods
        function this = DynLinTimoshenkoSystem(model)
            this = this@models.BaseDynSystem(model);
            
            this.MaxTimestep = [];
            
            %% Modellparameter
            % d1 = 0.3 D�mpfungsfaktor vor Massenmatrix (Luftwiderstand)
            this.addParam('d1',[0 2],10);
            % d2 = .01 D�mpfungsfaktor vor Steifigkeitsmatrix (Materiald�mpfung)
            this.addParam('d2',[0 .1],10);
            % Add the 3rd parameter, the local gravity factor
            this.addParam('gravity constant',[-20 20],10);
            
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
            
            data = this.Model.data;
            
            %% Assemblieren der Massenmatrix & Steifigkeitsmatrix f�r t=0
            % Mass matrix (M is actually the full M (including dirichlet
            % nodes), but is projected at the end of this method to the
            % size of actual DoFs)
            e = this.Model.Elements;
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
            M = sparse(i,j,M,7*nk,7*nk);
            % Project M to effectively needed entries
            M = M(this.Model.free,this.Model.free);
            this.M = dscomponents.ConstMassMatrix(blkdiag(speye(size(M)), M));
            this.M_small = M;
            % Create full K0 matrix
            this.K0 = sparse(i,j,K,7*nk,7*nk);
            
            dim = length(this.Model.free);
            
            %% Initial values
            x0 = zeros(2*dim,1);
            this.x0 = dscomponents.ConstInitialValue(x0);
                        
            %% Model input (gravity + const term)
            % Affine-parametric matrix models const. 
            % B(t,\mu) = \mu_3*[0 F] + 1*[f_const 0 0 0], u(t) = [1 g(t)], g\in\R^3
            B = dscomponents.AffLinInputConv;
            
            %% 1*[f_const 0 0 0]
            % Neumann forces (computed in Model.preprocess_data)
            f_const = this.Model.f_neum;
            f_const = f_const(this.Model.free);
            % Dirichlet forces only for linear model
            if ~this.Model.NonlinearModel
                % Reconstruct fake u vector which has entries only at
                % dirichlet points. Used for computation of constant forces due
                % to dirichlet values.
                u = zeros(7*nk,1);
                u(this.Model.dir_u) = this.Model.u_dir;
                f_dir = -this.K0*u;
                f_const = f_const + f_dir(this.Model.free);
            end
            Fc = sparse(2*dim,4);
            Fc(dim+1:end,4) = f_const;
            B.addMatrix('1',Fc);
            
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
            F = F(this.Model.free,:);
            F = [sparse(size(F,1),3); F];
            % Add zero column to front
            F = [sparse(size(F,1),1) F];
            B.addMatrix('mu(3)',F);
            B.CoeffClass = 'InputConvCoeffs';
            
            this.B = B;
            
            %% Output setup (dim = 3*space + 3*velo + 1*temp)
            %xdim = (dim - (dim/7))/2;
            nf = length(this.Model.free);
            C = speye(nf);
            % Remove velocity entries
            %C(repmat(logical([0 0 0 1 1 1]'),nf/6,1),:) = [];
            C = [C 0*C];
            this.C = dscomponents.LinearOutputConv(C);
            % Export settings
            je = this.Model.JKerMorExport;
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
        end
    end
end