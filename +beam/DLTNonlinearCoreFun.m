classdef DLTNonlinearCoreFun < models.beam.DLTBaseCoreFun
% DLTNonlinearCoreFun: 
%
%
%
% @author Daniel Wirtz @date 2011-12-07
%
% @new{0,6,dw,2011-12-07} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(Access=private)
        R_big;
        x_big = [];
        J_big;
    end
    
    methods
        function this = DLTNonlinearCoreFun(system)
            this = this@models.beam.DLTBaseCoreFun(system);
        end
        
        function fx = evaluateCoreFun(this, x, t, mu)%#ok
            % Daten nicht aktuell? (check inside)
            this.updatexJ(x);
            % Funktion f auswerten
            fx = this.f_big - this.B_big*x - this.R_big;
        end
        
        function J = getStateJacobian(this, x, t, mu)%#ok
            % Daten nicht aktuell? (check inside)
            this.updatexJ(x);
            % Jacobi-Matrix ausgeben
            J = this.J_big;
        end
        
        function pr = project(this, V, W)%#ok
            
        end
        
%         function prepareConstants(this, ~)
%             prepareConstants@models.beam.DLTBaseCoreFun(this);
%             
%             % Assemble constant part of B_big matrix (varying K for
%             % nonlinear case is added in evaluate)
%             
%         end
    end
    
    methods(Access=protected)
        function f_eff = getf_eff(~, f, ~, ~)
            f_eff = f;
        end
        
        function B = getB_big(~, ~, C)
            null = zeros(size(C));
            B = [null, -eye(size(C)); null, C];
        end
    end
    
    methods(Access=private)
        function updatexJ(this, x)
            m = this.sys.Model;
            % Updates the x, J and R big versions (if needed)
            if isempty(this.x_big) || (norm(this.x_big - x) > 1e-8)
                [R_F, K_F] = this.weak_form(m.data, m.RO, m.KR, m.FH, x);

                R = R_F(this.free);
                K = K_F(this.free, this.free);

                this.R_big = [zeros(size(R)); R];

                null = zeros(size(K));
                this.J_big = -1 * (this.B_big + [null, null; K, null]);
                this.x_big = x;
            end
        end
        
        function [R, K] = weak_form(this, data, RO, KR, FH, u)
            % Nichtlineare Funktion (Schwache Form, R-f_s) und ihre Ableitung (nach den Freiheitsgraden) (tangentielle Steifigkeitsmatrix, K)

            % Steifigkeitsmatrix für u und T
            K = sparse(7 * data.num_knots, 7 * data.num_knots);
            % Residuumsvektor
            R = sparse(7 * data.num_knots, 1);

            % Tangentiale Steifigkeitsmatrix aufstellen und schwache Form mit der aktuellen Verschiebung auswerten
            for i = 1:data.num_elem_RO

                index_0_glob = 7*data.knot_index(RO(i).p(1))-6 : 7*data.knot_index(RO(i).p(1))-1;
                index_L_glob = 7*data.knot_index(RO(i).p(2))-6 : 7*data.knot_index(RO(i).p(2))-1;
                index_glob = [index_0_glob index_L_glob];

                u_e = u(index_glob);

                [K_lok, R_lok] = this.beam_loc_matrix_tangential(u_e, RO(i).l, RO(i).T, RO(i).c);

                K(index_glob, index_glob) = K(index_glob, index_glob) + K_lok;
                R(index_glob) = R(index_glob) + R_lok;
            end

            for i = 1:data.num_elem_KR

                index_0_glob = 7*data.knot_index(KR(i).p(1))-6 : 7*data.knot_index(KR(i).p(1))-1;
                index_L_glob = 7*data.knot_index(KR(i).p(2))-6 : 7*data.knot_index(KR(i).p(2))-1;
                index_glob = [index_0_glob index_L_glob];

                u_e = u(index_glob);

                [K_lok, R_lok] = this.circle_loc_matrix_tangential(u_e, KR(i).R, KR(i).T_block1, KR(i).T_block2, KR(i).B, KR(i).angle, KR(i).c);

                K(index_glob, index_glob) = K(index_glob, index_glob) + K_lok;
                R(index_glob) = R(index_glob) + R_lok;
            end

            for i = 1:data.num_elem_FH

                index_0_glob = 7*data.knot_index(FH(i).p(1))-6 : 7*data.knot_index(FH(i).p(1))-4;
                index_L_glob = 7*data.knot_index(FH(i).p(2))-6 : 7*data.knot_index(FH(i).p(2))-4;
                index_glob = [index_0_glob index_L_glob];

                u_e = u(index_glob);

                [K_lok, R_lok] = this.truss_loc_matrix_tangential(u_e, FH(i).l, FH(i).T, FH(i).c);

                K(index_glob, index_glob) = K(index_glob, index_glob) + K_lok;
                R(index_glob) = R(index_glob) + R_lok;
            end
        end
        
        function [K, R, U_pot] = beam_loc_matrix_tangential(this, u, L, T_block, c)

            % Wertet tangentielle Steifigkeitsmatrix, "Residuumsvektor" (= Vektor der virtuellen internen Energie) und potenzielle Energie aus

            % c kodiert Stoffparameter:
            % c1            = E*I
            % c2 = c1^2     = (E*I)^2
            % c3            = G*As*L
            % c4 = c3^2     = (G*As*L)^2
            % c5            = G*As*L^2
            % c6 = c5^2     = (G*As*L^2)^2
            % c7 = c1*c5    = E*I*G*As*L^2
            % c8 = c1*G*As  = E*I*G*As
            % c9 = rho*A
            % c10= q = rho*A*Ortsfaktor
            % c11 = rho*I
            % c12 = rho*A*l/6; 
            % c13 = E*A/l;
            % c14 = G*I_t/l

            % Transformationsmatrix: natürliche Koords -> glob. Koords
            % Sortierung der Variablen: 
            % u0, v0, w0, phi0, psi0, theta0, u1, v1, w1, phi1, psi1, theta1
            %  1   2   3     4     5       6   7   8   9    10    11      12                        

            T = zeros(12);
            T([1 2 3], [1 2 3]) = T_block;         % Verschiebung Anfangsknoten
            T([7 8 9], [7 8 9]) = T_block;         % Verschiebung Endknoten
            T([4 5 6], [4 5 6]) = T_block;         % Winkel Anfangsknoten
            T([10 11 12], [10 11 12]) = T_block;   % Winkel Endknoten

            % Potenzielle Energie
            U_pot = 0;
            % Residuumsvektor
            R = zeros(12,1);
            % Tangentielle lokale Steifigkeitsmatrix
            K = zeros(12);

            % lokale Verschiebungen des Elements
            u_lok  = T' * u;

            % Stoff-Matrizen aufsetzen
            D = diag([c(13)*L, c(3)/L, c(3)/L]);
            E = diag([c(14)*L, c(1), c(1)]);


            %% Stützstellen und Gewichte für Gaußquadratur
            num_gauss_points = 3;
            switch (num_gauss_points)
                case 1
                    % Exakt bis Grad 1
                    s = 0.5*L;
                    w = L;
                case 2
                    % Exakt bis Grad 3
                    s = 0.5*L * [ (1-sqrt(1/3)), (1+sqrt(1/3)) ];
                    w = 0.5*L * [1 1];
                case 3
                    % Exakt bis Grad 5
                    s = 0.5*L * [ (1-sqrt(3/5)), 1, (1+sqrt(3/5)) ];
                    w = 0.5*L * 1/9 * [5 8 5];
                case 4
                    % Exakt bis Grad 7
                    s = 0.5*L * [ ( 1 - sqrt( 3/7 + 2/7*sqrt(6/5) ) ), 
                                  ( 1 - sqrt( 3/7 - 2/7*sqrt(6/5) ) ), 
                                  ( 1 + sqrt( 3/7 - 2/7*sqrt(6/5) ) ), 
                                  ( 1 + sqrt( 3/7 + 2/7*sqrt(6/5) ) )];
                    w = 0.5*L * 1/36 * [18-sqrt(30), 18+sqrt(30), 18+sqrt(30), 18-sqrt(30)];
                case 5
                    % Exakt bis Grad 9
                    s = 0.5*L * [ 1 - 1/3 * sqrt( 5 + 2*sqrt(10/7) ),
                                  1 - 1/3 * sqrt( 5 - 2*sqrt(10/7) ),
                                  1,
                                  1 + 1/3 * sqrt( 5 - 2*sqrt(10/7) ),
                                  1 + 1/3 * sqrt( 5 + 2*sqrt(10/7) )];
                    w = 0.5*L * 1/900 * [322 - 13*sqrt(70), 322 + 13*sqrt(70), 512, 322 + 13*sqrt(70), 322 - 13*sqrt(70)];
                otherwise
                    % Exakt bis Grad 5
                    s = 0.5*L * [ (1-sqrt(3/5)), 1, (1+sqrt(3/5)) ];
                    w = 0.5*L * 1/9 * [5 8 5];
            end


            % Kreuzproduktsmatrix e1 x v = S_e1 * v
            S_e1 = [0 0 0; 0 0 -1; 0 1 0];
            
            for i=1:length(s)
                % Auswerten der Ansatzfunktionen und deren Ableitungen
                N = models.beam.StraightBeam.beam_shape_functions(s(i), L, c);
                N_prime = models.beam.StraightBeam.beam_shape_functions_derivative(s(i), L, c);
            %     N1 = N(1:3,:);          % Verschiebungsansätze
                N2 = N(4:6,:);          % Verdrehungsansätze
                N1_prime = N_prime(1:3,:);
                N2_prime = N_prime(4:6,:);

                % Vorberechnungen zur Auswertung der Verzerrungen
                E_lin = N1_prime + S_e1 * N2;
                E_skw = N1_prime - S_e1 * N2;

                E_nl_1 = 0.125 * ( E_skw(2,:)' * E_skw(2,:) + E_skw(3,:)' * E_skw(3,:) );
                E_nl_2 = 0.5 * N2(1,:)' * E_skw(3,:);
                E_nl_3 = 0.5 * N2(1,:)' * E_skw(2,:);

                E_1 = 2 * E_nl_1;   % = E_nl_1 + E_nl_1';
                E_2 = E_nl_2 + E_nl_2';
                E_3 = E_nl_3 + E_nl_3';

                E_nonlin  = [u_lok' * E_nl_1; u_lok' * E_nl_2; u_lok' * E_nl_3];
                dE_nonlin = [u_lok' * E_1; u_lok' * E_2; u_lok' * E_3];
                dE = (E_lin + dE_nonlin);

                % Verzerrung
                eps = (E_lin + E_nonlin) * u_lok;

                % Krümmung
                kap = N2_prime * u_lok;

                % Kraft
                F = D * eps;
                % Moment
                M = E * kap;

                % Integrations-Update potenzielle Energie
                U_pot = U_pot + w(i) * ( eps' * F +  kap' * M);
                % Integrations-Update Schwache Form
                R = R + w(i) * ( dE' * F +  N2_prime' * M);
                % Integrations-Update Steifigkeitsmatrix
                K = K + w(i) * ( dE' * D * dE + N2_prime' * E * N2_prime + F(1)*E_1 + F(2)*E_2 + F(3)*E_3 );


            end

            % Transformation in globales Koordinatensystem (nicht nötig bei Skalaren)
            K = T * K * T';
            R = T * R;
        end
        
        function [K, R] = truss_loc_matrix_tangential(this, u, L, T_block, c)

            % Berechnet lokale tangentiale Steifigkeits- und Massenmatrix eines Stabes
            % c kodiert Stoffparameter:
            % c1 = E*A/L      (Federhärte)
            % c2 = rho*A*L/6

            % Transformationsmatrix: natürliche Koords -> glob. Koords
            % Sortierung der Variablen: 
            % u0, v0, w0, u1, v1, w1
            %  1   2   3   4   5   6                      

            T = zeros(6);
            T([1 2 3], [1 2 3]) = T_block;      % Verschiebung links
            T([4 5 6], [4 5 6]) = T_block;      % Verschiebung rechts

            % lokale Verschiebungen des Elements
            u_lok  = T' * u;
            % lokale Ableitungen
            u_lok_prime = (u_lok(4:6) - u_lok(1:3)) / L;
            % Lokales Spannungsmaß
            E_11 = u_lok_prime(1) + 0.5 * u_lok_prime' * u_lok_prime;
            S_x = L*c(1) * E_11;

            R = S_x * [-u_lok_prime - [1;0;0]; u_lok_prime + [1;0;0]];

            A1 = c(1) * (u_lok_prime + [1;0;0]) * (u_lok_prime + [1;0;0])';
            A2 = S_x/L * eye(3);
            B = A1 + A2;
            K = [B -B; -B B];

            K = T * K * T';
            R = T * R;
        end
        
        function [K, R, U_pot] = circle_loc_matrix_tangential(this, u, radius, T_block1, T_block2, B, angle, c)

            % Wertet tangentielle Steifigkeitsmatrix, "Residuumsvektor" (= Vektor der virtuellen internen Energie) und potenzielle Energie eines
            % Timoshenko-Balkens (Kreisbogen mit Öffnungswinkel [angle]) (numerisch mit speziell gewonnenen Ansatzfunktionen) aus

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


            % Potenzielle Energie
            U_pot = 0;
            % Residuumsvektor
            R = zeros(12,1);
            % Tangentielle lokale Steifigkeitsmatrix
            K = zeros(12);

            % lokale Verschiebungen des Elements
            u_lok  = T' * u;

            % Stoff-Matrizen aufsetzen
            D = blkdiag(c(3), c(4), c(4));
            E = blkdiag(c(2), c(1), c(1));

            L = radius * angle;
            %% Stützstellen und Gewichte für Gaußquadratur
            num_gauss_points = 3;
            switch (num_gauss_points)
                case 1
                    % Exakt bis Grad 1
                    s = 0.5*L;
                    w = L;
                case 2
                    % Exakt bis Grad 3
                    s = 0.5*L * [ (1-sqrt(1/3)), (1+sqrt(1/3)) ];
                    w = 0.5*L * [1 1];
                case 3
                    % Exakt bis Grad 5
                    s = 0.5*L * [ (1-sqrt(3/5)), 1, (1+sqrt(3/5)) ];
                    w = 0.5*L * 1/9 * [5 8 5];
                case 4
                    % Exakt bis Grad 7
                    s = 0.5*L * [ ( 1 - sqrt( 3/7 + 2/7*sqrt(6/5) ) ), 
                                  ( 1 - sqrt( 3/7 - 2/7*sqrt(6/5) ) ), 
                                  ( 1 + sqrt( 3/7 - 2/7*sqrt(6/5) ) ), 
                                  ( 1 + sqrt( 3/7 + 2/7*sqrt(6/5) ) )];
                    w = 0.5*L * 1/36 * [18-sqrt(30), 18+sqrt(30), 18+sqrt(30), 18-sqrt(30)];
                case 5
                    % Exakt bis Grad 9
                    s = 0.5*L * [ 1 - 1/3 * sqrt( 5 + 2*sqrt(10/7) ),
                                  1 - 1/3 * sqrt( 5 - 2*sqrt(10/7) ),
                                  1,
                                  1 + 1/3 * sqrt( 5 - 2*sqrt(10/7) ),
                                  1 + 1/3 * sqrt( 5 + 2*sqrt(10/7) )];
                    w = 0.5*L * 1/900 * [322 - 13*sqrt(70), 322 + 13*sqrt(70), 512, 322 + 13*sqrt(70), 322 - 13*sqrt(70)];
                otherwise
                    % Exakt bis Grad 5
                    s = 0.5*L * [ (1-sqrt(3/5)), 1, (1+sqrt(3/5)) ];
                    w = 0.5*L * 1/9 * [5 8 5];
            end

            %%
            % Kreuzproduktsmatrix e1 x v = S_e1 * v
            S_e1 = [0 0 0; 0 0 -1; 0 1 0];
            % Matrizen, die die lokale Verschiebung auf ihre (über den Kreisbogen
            % konstante!) *LINEARE* Verzerrung und Krümmung abbilden
            B3 = B(7:9,:);
            B4 = B(10:12,:);

            for i=1:length(s)
                % Auswerten der Ansatzfunktionen
                N = models.beam.CurvedBeam.circle_shape_functions(radius, s(i), B);
            %     N1 = N(1:3,:);          % Verschiebungsansätze
                N2 = N(4:6,:);          % Verdrehungsansätze

                % Vorberechnungen zur Auswertung der Verzerrungen
                E_lin = B3;      % = (für gerade Balken) N1_prime + S_e1 * N2;
                E_skw = E_lin - 2 * S_e1 * N2;  % = N1_prime - S_e1 * N2;

                E_nl_1 = 0.125 * ( E_skw(2,:)' * E_skw(2,:) + E_skw(3,:)' * E_skw(3,:) );
                E_nl_2 = 0.5 * N2(1,:)' * E_skw(3,:);
                E_nl_3 = 0.5 * N2(1,:)' * E_skw(2,:);

                E_1 = 2 * E_nl_1;   % = E_nl_1 + E_nl_1';
                E_2 = E_nl_2 + E_nl_2';
                E_3 = E_nl_3 + E_nl_3';

                E_nonlin  = [u_lok' * E_nl_1; u_lok' * E_nl_2; u_lok' * E_nl_3];
                dE_nonlin = [u_lok' * E_1; u_lok' * E_2; u_lok' * E_3];
                dE = (E_lin + dE_nonlin);

                % Verzerrung
                eps = (E_lin + E_nonlin) * u_lok;    
                % Krümmung
                kap = B4 * u_lok;       % = N2_prime * u_lok;

                % Kraft
                F = D * eps;
                % Moment
                M = E * kap;

                % Integrations-Update potenzielle Energie
                U_pot = U_pot + w(i) * ( eps' * F +  kap' * M);
                % Integration
                R = R + w(i) * ( dE' * F +  B4' * M);
                K = K + w(i) * ( dE' * D * dE + B4' * E * B4 + F(1)*E_1 + F(2)*E_2 + F(3)*E_3 );

            end

            % Transformation in globales Koordinatensystem (nicht nötig bei Skalaren)
            K = T * K * T';
            R = T * R;
        end
    end
    
end