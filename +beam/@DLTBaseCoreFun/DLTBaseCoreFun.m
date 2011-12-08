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
%                 [M_lok, K_lok] = this.loc_matrix_beam(RO(i).l, RO(i).T, RO(i).c);
%                 f_lok = this.loc_matrix_beam_force(RO(i).l, RO(i).T, RO(i).q_lok);
                M_lok = RO(i).getLocalMassMatrix;
                K_lok = RO(i).getLocalStiffnessMatrix;
                f_lok = RO(i).getLocalForce(m.Gravity);
                index_0_lok = [1 5 9 3 10 6];
                index_l_lok = [2 7 11 4 12 8];

                p = RO(i).PointsIdx;
                index_0_glob = 7*data.knot_index(p(1))-6 : 7*data.knot_index(p(1))-1;
                index_l_glob = 7*data.knot_index(p(2))-6 : 7*data.knot_index(p(2))-1;

                index_glob = [index_0_glob index_l_glob];
                index_lok = [index_0_lok index_l_lok];

                M(index_glob, index_glob) = M(index_glob, index_glob) + M_lok(index_lok, index_lok); %#ok<*SPRIX>
                K(index_glob, index_glob) = K(index_glob, index_glob) + K_lok(index_lok, index_lok);
                f_const(index_glob) = f_const(index_glob) + f_lok(index_lok);
            end
            
            KR = m.KR;
            for i = 1:data.num_elem_KR
                [M_lok, K_lok] = this.loc_matrix_circle_special(KR(i).R, KR(i).T_block1, KR(i).T_block2, KR(i).B, KR(i).angle, KR(i).c);
                f_lok = this.loc_matrix_circle_special_force(KR(i).R, KR(i).Fren, KR(i).T_block1, KR(i).T_block2, KR(i).B, KR(i).angle, KR(i).q_lok);

                index_0_glob = 7*data.knot_index(KR(i).p(1))-6 : 7*data.knot_index(KR(i).p(1))-1;
                index_l_glob = 7*data.knot_index(KR(i).p(2))-6 : 7*data.knot_index(KR(i).p(2))-1;

                index_glob = [index_0_glob index_l_glob];

                M(index_glob, index_glob) = M(index_glob, index_glob) + M_lok;
                K(index_glob, index_glob) = K(index_glob, index_glob) + K_lok;
                f_const(index_glob) = f_const(index_glob) + f_lok;

            end

            FH = m.FH;
            for i = 1:data.num_elem_FH

                [M_lok, K_lok] = this.loc_matrix_truss(FH(i).T, FH(i).c);
                f_lok = this.loc_matrix_truss_force(FH(i).l, FH(i).T, FH(i).q_lok);

                index_0_glob = 7*data.knot_index(FH(i).p(1))-6 : 7*data.knot_index(FH(i).p(1))-4;
                index_l_glob = 7*data.knot_index(FH(i).p(2))-6 : 7*data.knot_index(FH(i).p(2))-4;
                index_glob = [index_0_glob index_l_glob];

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
            
%             this.free_uv = [free_u free_u + 7 * data.num_knots];
%             this.free_u = free_u;
        end
    end
    
    methods(Access=protected)
        f_eff = getf_eff(this, f, K, u);
        
        B = getB_big(this, K, C);
    end
    
    methods(Access=private)
        function [M, K] = loc_matrix_circle_special(this, R, T_block1, T_block2, B, angle, c)

            % Berechnet lokale Steifigkeits- und Massenmatrix eines gekrümmten
            % Timoshenko-Balkens (Viertelkreis) (numerisch mit speziell gewonnenen Ansatzfunktionen)

            % c1 = E*I
            % c2 = G*It
            % c3 = E*A
            % c4 = G*As
            % c5 = rho*A
            % c6 = rho*I
            % c7 = rho*It

            % Transformationsmatrix für lokale in globale Größen
            % T_block = [e_x e_y e_z];
            % Transformationsmatrix für Größen am Anfangspunkt (Frenet_Basis in globale
            % Koordinaten umgerechnet
            % T_block1 = T_block * T_Fren(0);
            % Transformationsmatrix für Größen am Endpunkt
            % T_block2 = T_block * T_Fren(angle);

            % Transformationsmatrix: natürliche Koords -> glob. Koords
            % Sortierung der Variablen: 
            % u0, v0, w0, phi0, psi0, theta0, u1, v1, w1, phi1, psi1, theta1
            %  1   2   3     4     5       6   7   8   9    10    11      12                        

            T = zeros(12);
            T([1 2 3], [1 2 3]) = T_block1;         % Verschiebung Anfangsknoten
            T([7 8 9], [7 8 9]) = T_block2;         % Verschiebung Endknoten
            T([4 5 6], [4 5 6]) = T_block1;         % Winkel Anfangsknoten
            T([10 11 12], [10 11 12]) = T_block2;   % Winkel Endknoten


            % Ansatzfunktionen N und Verbindungsmatrix für einen Viertelkreis holen

            % Matrizen aufsetzen
            D = blkdiag(c(3), c(4), c(4));
            E = blkdiag(c(2), c(1), c(1));
            rhoJ = blkdiag(c(7), c(6), c(6));

            M = zeros(12);

            L = R*angle;
            % Stützstellen für Gaußquadratur
            s = 0.5*L * [ (1-sqrt(3/5)), 1, (1+sqrt(3/5)) ];
            % Gewichte für Gaußquadratur
            w = 0.5*L * 1/9 * [5 8 5];

            for i=1:3
                N = models.beam.CurvedBeam.circle_shape_functions(R, s(i), B);
                N1 = N(1:3,:);
                N2 = N(4:6,:);    
                M = M + w(i) * ( c(5)*N1'*N1 + N2'*rhoJ*N2 );
            end

            B3 = B(7:9,:);
            B4 = B(10:12,:);

            K = L * (B4' * E * B4 + B3' * D * B3);

            M = T * M * T';
            K = T * K * T';
        end
        
        function f = loc_matrix_circle_special_force(this, R, T_Fren, T_block1, T_block2, B, angle, q_lok)

            % Berechnet lokale Steifigkeits- und Massenmatrix eines gekrümmten
            % Timoshenko-Balkens (Viertelkreis) (numerisch mit speziell gewonnenen Ansatzfunktionen)

            % c1 = E*I
            % c2 = G*It
            % c3 = E*A
            % c4 = G*As
            % c5 = rho*A
            % c6 = rho*I
            % c7 = rho*It

            % Transformationsmatrix für lokale in globale Größen
            % T_block = [e_x e_y e_z];
            % Transformationsmatrix für Größen am Anfangspunkt (Frenet_Basis in globale
            % Koordinaten umgerechnet
            % T_block1 = T_block * T_Fren(0);
            % Transformationsmatrix für Größen am Endpunkt
            % T_block2 = T_block * T_Fren(angle);

            % Transformationsmatrix: natürliche Koords -> glob. Koords
            % Sortierung der Variablen: 
            % u0, v0, w0, phi0, psi0, theta0, u1, v1, w1, phi1, psi1, theta1
            %  1   2   3     4     5       6   7   8   9    10    11      12                        

            T = zeros(12);
            T([1 2 3], [1 2 3]) = T_block1;         % Verschiebung Anfangsknoten
            T([7 8 9], [7 8 9]) = T_block2;         % Verschiebung Endknoten
            T([4 5 6], [4 5 6]) = T_block1;         % Winkel Anfangsknoten
            T([10 11 12], [10 11 12]) = T_block2;   % Winkel Endknoten

            f = zeros(12,1);

            L = R*angle;
            % Stützstellen für Gaußquadratur
            s = 0.5*L * [ (1-sqrt(3/5)), 1, (1+sqrt(3/5)) ];
            % Gewichte für Gaußquadratur
            w = 0.5*L * 1/9 * [5 8 5];
            for i=1:3
                N = models.beam.CurvedBeam.circle_shape_functions(R, s(i), B);
                f = f + w(i) * N(1:3,:)'*(T_Fren(s(i)/R)'*q_lok);
            end

            f = T * f;
        end
        
        function [M, K] = loc_matrix_truss(this, T_block, c)
            % Berechnet lokale Steifigkeits- und Massenmatrix eines Stabes
            % c kodiert Stoffparameter:
            % c1 = E*A/L      (Federhärte)
            % c2 = rho*A*L/6

            % Transformationsmatrix: natürliche Koords -> glob. Koords
            % Sortierung der Variablen: 
            % u0, v0, w0, u1, v1, w1
            %  1   2   3   4   5   6                      

            T = zeros(6);
            T([1 2 3], [1 2 3]) = T_block;      % Verschiebung links
            T([4 5 6], [4 5 6]) = T_block;    % Verschiebung rechts

            K = zeros(6);
            K([1 4], [1 4]) = c(1) * [1 -1; -1 1];

            M = zeros(6);
            M([1 4], [1 4]) = c(2) * [2 1; 1 2];

            M = T * M * T';
            K = T * K * T';
        end
        
        function f = loc_matrix_truss_force(this, l, T_block, q_lok)
            % Berechnet lokale Steifigkeits- und Massenmatrix eines Stabes
            % c kodiert Stoffparameter:
            % c1 = E*A/L      (Federhärte)
            % c2 = rho*A*L/6

            % Transformationsmatrix: natürliche Koords -> glob. Koords
            % Sortierung der Variablen: 
            % u0, v0, w0, u1, v1, w1
            %  1   2   3   4   5   6                      

            T = zeros(6);
            T([1 2 3], [1 2 3]) = T_block;      % Verschiebung links
            T([4 5 6], [4 5 6]) = T_block;    % Verschiebung rechts

            f = zeros(6, 1);
            f([1 4]) = 0.5 * q_lok(1) * l * [1; 1];

            f = T * f;
        end
    end
end