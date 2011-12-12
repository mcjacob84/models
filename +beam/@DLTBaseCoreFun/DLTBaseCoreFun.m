classdef DLTBaseCoreFun < dscomponents.ACoreFun & dscomponents.IJacobian
% DynLinTimoshenkoCoreFun: 
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
    
    properties(SetAccess=private, GetAccess=protected)
        sys;
        % The stiffness matrix for t=0, used in the linear model and also
        % in the nonlinear model to compute C.
        K0;
    end
    
    properties(Access=protected)
        f_big;
    end
    
    methods
        function this = DLTBaseCoreFun(sys)
            this = this@dscomponents.ACoreFun;
            this.sys = sys;
            this.initialize;
        end
       
%         function copy = clone(this)%#ok
%             
%         end

        function initialize(this)
            % do constant preps, i.e. assign A
            % Massenmatrix (u und T)
            m = this.sys.Model;
            data = m.data;
            
            %% Assemblieren der globalen matrizen
%             % Steifigkeitsmatrix für u und T
%             K = sparse(7 * data.num_knots, 7 * data.num_knots);
%             % Kraftvektor und Wärmequellen 
%             f_const = sparse(7 * data.num_knots, 1);
            e = this.sys.Model.Elements;
            i = []; j = [];
            K = []; 
            fi = []; f = [];
            for k = 1:length(e)
                K_lok = e{k}.getLocalStiffnessMatrix;
                f_lok = e{k}.getLocalForce(m.Gravity);
                index_glob = e{k}.getGlobalIndices;

                [li,lj] = meshgrid(index_glob);
                i = [i; li(:)];%#ok<*AGROW>
                j = [j; lj(:)];
                K = [K; K_lok(:)];
                fi = [fi; index_glob'];%#ok<*AGROW>
                f = [f; f_lok(:)];
%                 K(index_glob, index_glob) = K(index_glob, index_glob) + K_lok;
%                 f_const(index_glob) = f_const(index_glob) + f_lok;
            end
            this.K0 = sparse(i,j,K,7*data.num_knots,7*data.num_knots);
            f_const = sparse(fi,ones(size(fi)),f,7*data.num_knots,1);
            
            % Add Neumann forces (computed in Model.preprocess_data)
            f_const = f_const + m.f_neum;
            
            % Dirichlet forces will be added directly in LinearFun.initialize, but
            % as there are none in NonlinearFun there is nothing to do here
            
            f_const = f_const(m.free);
            this.f_big = [zeros(size(f_const)); f_const];
            
            %% Jacobian matrix sparsity pattern
            % Create a fake big system matrix here to obtain the sparsity
            % pattern.
            s = length(m.free);
            K = this.K0(m.free,m.free);
            null = sparse(s,s);
            [i,j] = find([null, -speye(s,s); K, K]);
            this.JSparsityPattern = sparse(i,j,ones(size(i)),2*s,2*s);
        end
    end
end