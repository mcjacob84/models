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
    end
        
    methods
        function this = DynLinTimoshenkoSystem(model)
            this = this@models.BaseDynSystem(model);
            
            this.MaxTimestep = [];
            
            %% Modellparameter
            % d1 = 0.3 Dämpfungsfaktor vor Massenmatrix (Luftwiderstand)
            this.addParam('d1',[0 2],5);
            % d2 = .01 Dämpfungsfaktor vor Steifigkeitsmatrix (Materialdämpfung)
            this.addParam('d2',[0 2],5);
            
            %% Assemblieren der Massenmatrix
            % Mass matrix (M is actually the full M (including dirichlet
            % nodes), but is projected at the end of this method to the
            % size of actual DoFs)
%             M2 = sparse(7 * model.data.num_knots, 7 * model.data.num_knots);
            e = this.Model.Elements;
            i = []; j = [];
            M = []; 
            for k = 1:length(e)
                M_lok = e{k}.getLocalMassMatrix;
                index_glob = e{k}.getGlobalIndices;

                [li,lj] = meshgrid(index_glob);
                i = [i; li(:)]; %#ok<*AGROW>
                j = [j; lj(:)];
                M = [M; M_lok(:)];
                
%                 M2(index_glob, index_glob) = M2(index_glob, index_glob) + M_lok;
            end
            nk = model.data.num_knots;
            M = sparse(i,j,M,7*nk,7*nk);
            %% Mass matrix for system
            % Project M to effectively needed entries
            M = M(model.free,model.free);
            this.M = dscomponents.ConstMassMatrix(blkdiag(eye(size(M)), M));
            this.M_small = M;
            
            %% Initial values
            x0 = zeros(2*length(model.free),1);
            this.x0 = dscomponents.ConstInitialValue(x0);
            
            %% Input (gravity)
            % Add the 3rd parameter, the local gravity factor
            this.addParam('gravity constant',[0 20],1);

            i = []; j = [];
            F = []; 
            for k = 1:length(e)
                F_lok = e{k}.getLocalForceMatrix;
                index_glob = e{k}.getGlobalIndices;

                [li,lj] = meshgrid(index_glob,1:3);
                i = [i; li(:)]; %#ok<*AGROW>
                j = [j; lj(:)];
                F = [F; F_lok(:)];
            end
            B = dscomponents.AffLinInputConv;
            F = sparse(i,j,F,7*nk,3);
            F = F(model.free,:);
            B.addMatrix('mu(3)',[sparse(size(F,1),3); F]);
            this.B = B;
            
            % Fake constant gravity
            this.Inputs{1} = @(t)[0; 0; -1];
            this.Inputs{2} = @(t)[-1; 0; 0];
            this.Inputs{3} = @(t)[0; -1; 0];
            % More "advanced" gravity
            this.Inputs{4} = @(t)[sin(t); 0; 0];
            T = model.T;
            this.Inputs{5} = @(t)[sin((t/T)*2*pi); cos((t/T)*2*pi); 0];
            this.Inputs{6} = @(t)[0; sin((t/T)*2*pi); cos((t/T)*2*pi)];
            
%             %% Output setup (dim = 3*space + 3*velo + 1*temp)
%             %xdim = (dim - (dim/7))/2;
%             d = sparse(1:dim,1:dim,1);
%             d(repmat(logical([0 0 0 1 1 1 0]'),dim/7,1),:) = [];
%             %d = reshape(repmat([1 1 1 0 0 0 1]', dim/7,1),[],1);
%             %d = [d; 0*d];
%             this.C = dscomponents.LinearOutputConv(d);
        end
    end
end