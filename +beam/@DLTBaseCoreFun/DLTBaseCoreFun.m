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
    
    properties(SetAccess=private)
        dir_u;
        dir_T;
        
        % The indices in the global state space vector of all points
        % including dirichlet points (dim=7: 3location, 3velocity & heat) 
        free;
        
        % Extracted Dirichlet values full u vector (using dir_u)
        u_dir;
    end

    properties(SetAccess=private, GetAccess=protected)
        B_big;
        f_big;
        sys;
    end
    
    methods
        function this = DLTBaseCoreFun(sys)
            this = this@dscomponents.ACoreFun;
            this.sys = sys;
        end
       
%         function copy = clone(this)%#ok
%             
%         end
        
        function prepareConstants(this, ~)
            % do constant preps, i.e. assign A
            % Massenmatrix (u und T)
            m = this.sys.Model;
            data = m.data;
            RO = m.RO;
            withHeat = false;
            
            % Mass matrix (M is actually the full M (including dirichlet
            % nodes), but is projected at the end of this method to the
            % size of actual DoFs)
            M = sparse(7 * data.num_knots, 7 * data.num_knots);
            % Steifigkeitsmatrix für u und T
            K = sparse(7 * data.num_knots, 7 * data.num_knots);
            % Kraftvektor und Wärmequellen 
            f_const = sparse(7 * data.num_knots, 1); 

            for i = 1:data.num_elem_RO
                M_lok = RO(i).getLocalMassMatrix;
                K_lok = RO(i).getLocalStiffnessMatrix;
                f_lok = RO(i).getLocalForce(m.Gravity);
                index_0_lok = [1 5 9 3 10 6];
                index_l_lok = [2 7 11 4 12 8];

%                 p = RO(i).PointsIdx;
%                 index_0_glob = 7*data.knot_index(p(1))-6 : 7*data.knot_index(p(1))-1;
%                 index_l_glob = 7*data.knot_index(p(2))-6 : 7*data.knot_index(p(2))-1;
% 
%                 index_glob = [index_0_glob index_l_glob];
                index_glob = RO(i).getGlobalIndices;
                index_lok = [index_0_lok index_l_lok];

                M(index_glob, index_glob) = M(index_glob, index_glob) + M_lok(index_lok, index_lok); %#ok<*SPRIX>
                K(index_glob, index_glob) = K(index_glob, index_glob) + K_lok(index_lok, index_lok);
                f_const(index_glob) = f_const(index_glob) + f_lok(index_lok);
            end
            
            KR = m.KR;
            for i = 1:data.num_elem_KR
                M_lok = KR(i).getLocalMassMatrix;
                K_lok = KR(i).getLocalStiffnessMatrix;
                f_lok = KR(i).getLocalForce(m.Gravity);
                
                %[M_lok, K_lok] = this.loc_matrix_circle_special(KR(i).R, KR(i).T_block1, KR(i).T_block2, KR(i).B, KR(i).angle, KR(i).c);
                %f_lok = this.loc_matrix_circle_special_force(KR(i).R, KR(i).Fren, KR(i).T_block1, KR(i).T_block2, KR(i).B, KR(i).angle, KR(i).q_lok);
%                 p = KR(i).PointsIdx;
%                 index_0_glob = 7*data.knot_index(p(1))-6 : 7*data.knot_index(p(1))-1;
%                 index_l_glob = 7*data.knot_index(p(2))-6 : 7*data.knot_index(p(2))-1;
% 
%                 index_glob = [index_0_glob index_l_glob];
                index_glob = KR(i).getGlobalIndices;
                M(index_glob, index_glob) = M(index_glob, index_glob) + M_lok;
                K(index_glob, index_glob) = K(index_glob, index_glob) + K_lok;
                f_const(index_glob) = f_const(index_glob) + f_lok;
            end

            FH = m.FH;
            for i = 1:data.num_elem_FH
                M_lok = FH(i).getLocalMassMatrix;
                K_lok = FH(i).getLocalStiffnessMatrix;
                f_lok = FH(i).getLocalForce(m.Gravity);

%                 index_0_glob = 7*data.knot_index(FH(i).PointsIdx(1))-6 : 7*data.knot_index(FH(i).PointsIdx(1))-4;
%                 index_l_glob = 7*data.knot_index(FH(i).PointsIdx(2))-6 : 7*data.knot_index(FH(i).PointsIdx(2))-4;
%                 index_glob = [index_0_glob index_l_glob];
                index_glob = FH(i).getGlobalIndices;
                
                M(index_glob, index_glob) = M(index_glob, index_glob) + M_lok;
                K(index_glob, index_glob) = K(index_glob, index_glob) + K_lok;
                f_const(index_glob) = f_const(index_glob) + f_lok;

%                 % Beiträge der Heat-equation assemblieren
%                 [M_lok, K_lok, f_lok] = loc_matrix_heat(FH(i).l, [0; 0]);
%                 index_glob = [7*data.knot_index(FH(i).p(1)) 7*data.knot_index(FH(i).p(2))];
%                 if (dynamic)
%                     M(index_glob, index_glob) = M(index_glob, index_glob) + FH(i).c_theta * M_lok;
%                 end
%                 K(index_glob, index_glob) = K(index_glob, index_glob) + FH(i).kappa * K_lok;    
            end

            %% Einbau der Randbedingungen
            % Balken
            for i = 1:length(data.neumann)
                if (data.knot_index(data.neumann(i)) == 0)
                    continue;
                end
                index_start = 7*data.knot_index(data.neumann(i))-6;
                f_const(index_start : index_start+6) = f_const(index_start : index_start+6) + data.neu_data(i, :)';
            end

            nodes = 1: 7*data.num_knots;
            nodes_T = 7 * (1:data.num_knots);
            nodes_u = setdiff(nodes, nodes_T);

            % Aufbau des vollen Dirichlet-Vektors
            % Vektor mit Indizes der Freiheitsgrade von u aufbauen
            u = zeros(7 * data.num_knots, 1);
            free = 1:7*data.num_knots;

            % Balken
            for i = 1:length(data.dirichlet)
                % Lokale Liste der Komponenten dieses Knotens, die Dirichlet sind
                dir_knoten_lok = (1:7) .* data.dir_data(i,8:14);
                dir_knoten_lok = setdiff(dir_knoten_lok, 0);

                % Wenn der Knoten in der Struktur nicht verwendet wird, auslassen
                if (data.knot_index(data.dirichlet(i)) == 0)
                    continue;
                end

                % Globale Liste der Komponenten dieses Knotens, die Dirichlet sind
                dir_knoten = 7*data.knot_index(data.dirichlet(i))-6 : 7*data.knot_index(data.dirichlet(i));
                dir_knoten = dir_knoten .* data.dir_data(i,8:14);
                dir_knoten = setdiff(dir_knoten, 0);

                % Dirichlet-Werte in Lösungsvektor eintragen und Liste der Freiheitsgrade verkleinern
                u(dir_knoten) = data.dir_data(i, dir_knoten_lok)';
                free = setdiff(free, dir_knoten);
            end

            %% Vektor der Freiheitsgrade aufteilen in Freiheitsgrade für Balken und Temperatur
            % Vektor für Balken (Temperatureinträge entfernen)
            if withHeat
                free_T = setdiff(free, free_u);
                free_u = setdiff(free, nodes_T);
            else
                free_u = setdiff(free,nodes_T);
                free = free_u;
            end
            this.free = free;
            
            this.dir_u = setdiff(nodes_u, free_u);
            if withHeat
                this.dir_T = setdiff(nodes_T, free_T);
            end
            % Store dirichlet values
            this.u_dir = u(this.dir_u);
            
            %% Mass matrix for system
            % Project M to effectively needed entries
            M = M(free,free);
            this.sys.M = dscomponents.ConstMassMatrix(blkdiag(eye(size(M)), M));
            
            % Call template method, actual computation in subclasses
            f_const_eff = this.getf_eff(f_const,K,u);
            
            f_const_eff = f_const_eff(free);
            this.f_big = [zeros(size(f_const_eff)); f_const_eff];
            
            % Project K now
            K = K(free,free);
            
            % Dämpfungsmodell 1: M a_t + (d1*M + d2*K) v_t + K u_t = f
            d1 = 0.3;   % Dämpfungsfaktor vor Massenmatrix (Luftwiderstand)
            d2 = .01;   % Dämpfungsfaktor vor Steifigkeitsmatrix (Materialdämpfung)
            C = (d1*M + d2*K); %#ok<*PROP>
            
            % Store K and C values so that the subclass prepareConstants
            % method can use them for assemblation of B_big matrices
            this.B_big = this.getB_big(K,C);
            
            %% Jacobian matrix sparsity pattern
            % Create a fake B_big here as the nonlinear getB_big returns
            % only the constant part of B_big. However, this does not
            % affect the sparsity pattern, so create it here.
            [i,j] = find([zeros(size(K)), -eye(size(K)); K, C]);
            this.JSparsityPattern = sparse(i,j,ones(size(i)));
        end
    end
    
    methods(Access=protected)
        f_eff = getf_eff(this, f, K, u);
        
        B = getB_big(this, K, C);
    end
end