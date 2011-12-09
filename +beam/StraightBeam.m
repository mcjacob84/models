classdef StraightBeam < models.beam.Beam
% StraightBeam: 
%
%
%
% @author Daniel Wirtz @date 2011-12-05
%
% @new{0,6,dw,2011-12-05} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods
        
        function this = StraightBeam(model, material, pointsidx)
            this = this@models.beam.Beam(model, material, pointsidx);
            this.initialize;
        end
        
        function M = getLocalMassMatrix(this)
            % Berechnet lokale Steifigkeits- und Massenmatrix eines Timoshenko-Balkens
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
            
            c = this.c;
            T_block = this.T;
            l = this.Length;

            % Transformationsmatrix: natürliche Koords -> glob. Koords
            % Schubwinkel beta wird nicht transformiert (eye(4))
            % Sortierung der Variablen: 
            % u0, u1, phi0, phi1, v0, theta0, v1, theta1, w0, psi0, w1, psi1
            %  1   2     3     4   5       6   7       8   9    10  11    12 
            T = zeros(12);
            T([1 5 9], [1 5 9]) = T_block;      % Verschiebung links
            T([2 7 11], [2 7 11]) = T_block;    % Verschiebung rechts
            T([3 10 6], [3 10 6]) = T_block;    % Winkel links
            T([4 12 8], [4 12 8]) = T_block;    % Winkel rechts

            M_u_block = c(12) * [2 1; 1 2];
            M_phi_block = c(15) * [2 1; 1 2];

            d1 = 12*(1680*c(2) + 294*c(7) + 13*c(6));
            d2 = 4*l^2*(126*c(2) + 21*c(7) + c(6));
            n1 = 2*l*(1260*c(2) + 231*c(7) + 11*c(6));
            n2 = l*(2520*c(2) + 378*c(7) + 13*c(6));
            b1 = 18*(560*c(2) + 84*c(7) + 3*c(6));
            b2 = 3*l^2*(168*c(2) + 28*c(7) + c(6));

            factor1 = c(9)*l / (420 * (12*c(1) + c(5))^2 );

            M_v =factor1 * [
                 d1  n1  b1 -n2;
                 n1  d2  n2  -b2;
                 b1  n2  d1 -n1;
                -n2  -b2 -n1  d2];

            M_w = factor1 * [
                 d1 -n1  b1  n2;
                -n1  d2 -n2  -b2;
                 b1 -n2  d1  n1;
                 n2  -b2  n1  d2];

            d1 = 36*c(4);
            d2 = 4*(360*c(2) + 15*c(7) + c(6));
            n1 = 3*c(3)*(-60*c(1) + c(5));
            b2 = 720*c(2) - 60*c(7) - c(6);
            factor2 = c(11)*l / (30 * (12*c(1) + c(5))^2 );

            M_theta = factor2 * [
                 d1  n1 -d1  n1;
                 n1  d2 -n1  b2;
                -d1 -n1  d1 -n1;
                 n1  b2 -n1  d2];
            M_psi = factor2 * [
                 d1 -n1 -d1 -n1;
                -n1  d2  n1  b2;
                -d1  n1  d1  n1;
                -n1  b2  n1  d2];

            M = [M_u_block      zeros(2,2)      zeros(2,4)  zeros(2,4);
                 zeros(2,2)     M_phi_block     zeros(2,4)  zeros(2,4);
                 zeros(4,2)     zeros(4,2)      M_v+M_theta zeros(4,4);
                 zeros(4,2)     zeros(4,2)      zeros(4,4)  M_w+M_psi];

            M = T * M * T';
        end
        
        function K = getLocalStiffnessMatrix(this)
            % Berechnet lokale Steifigkeitsmatrix eines Timoshenko-Balkens
            %
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
            
            c = this.c;
            T_block = this.T;
            l = this.Length;

            % Transformationsmatrix: natürliche Koords -> glob. Koords
            % Schubwinkel beta wird nicht transformiert (eye(4))
            % Sortierung der Variablen: 
            % u0, u1, phi0, phi1, v0, theta0, v1, theta1, w0, psi0, w1, psi1
            %  1   2     3     4   5       6   7       8   9    10  11    12 
            T = zeros(12);
            T([1 5 9], [1 5 9]) = T_block;      % Verschiebung links
            T([2 7 11], [2 7 11]) = T_block;    % Verschiebung rechts
            T([3 10 6], [3 10 6]) = T_block;    % Winkel links
            T([4 12 8], [4 12 8]) = T_block;    % Winkel rechts

            K_u_block = c(13) * [1 -1; -1 1];   % * E A/L
            K_phi_block = c(14)* [1 -1; -1 1];  % * G I_t/L

            d1 = 12*c(8);
            d2 = 4*c(1)*(3*c(1) + c(5));
            n1 = 6*l*c(8);
            b2 = 2*c(1)*(c(5) - 6*c(1));
            factor3 = 1 / ( l * (12*c(1) + c(5)) );

            K_v_block = factor3 * [
                d1  n1 -d1  n1;
                n1  d2 -n1  b2;
               -d1 -n1  d1 -n1;
                n1  b2 -n1  d2];

            K_w_block = factor3 * [
                d1 -n1 -d1 -n1;
               -n1  d2  n1  b2;
               -d1  n1  d1  n1;
               -n1  b2  n1  d2];

            K = [K_u_block      zeros(2,2)      zeros(2,4)  zeros(2,4);
                 zeros(2,2)     K_phi_block     zeros(2,4)  zeros(2,4);
                 zeros(4,2)     zeros(4,2)      K_v_block   zeros(4,4);
                 zeros(4,2)     zeros(4,2)      zeros(4,4)  K_w_block];

            K = T * K * T';
        end
        
        function f = getLocalForce(this, gravity)
            l = this.Length;
            q_lok = (this.c(9) + this.ROHR_q_plus) * (this.T' * gravity);
            
            f_u = q_lok(1) * 0.5 * l * [1; 1; 0; 0];
            f_v = q_lok(2) * l/12 * [6; l; 6; -l];
            f_w = q_lok(3) * l/12 * [6; -l; 6; l];
            
            T_block = this.T;
            T = zeros(12);
            T([1 5 9], [1 5 9]) = T_block;      % Verschiebung links
            T([2 7 11], [2 7 11]) = T_block;    % Verschiebung rechts
            T([3 10 6], [3 10 6]) = T_block;    % Winkel links
            T([4 12 8], [4 12 8]) = T_block;    % Winkel rechts
            
            f = T * [f_u; f_v; f_w];
        end
        
        function [K, R, U_pot] = getLocalTangentials(this, u)
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
            
            L = this.Length;
            T_block = this.T;
            c = this.c;

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
                N_prime = this.beam_shape_functions_derivative(s(i));
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
        
        function N = beam_shape_functions_derivative(this, s)
            % Wertet (analytische) Basisfunktion für einen geraden Balken der Länge L und den Stoffkonstanten c an der Stelle s aus.
            % c1            = E*I
            % c2 = c1^2     = (E*I)^2
            % c3            = G*As*L
            % c4 = c3^2     = (G*As*L)^2
            % c5 = c3*L     = G*As*L^2
            % c6 = c5^2     = (G*As*L^2)^2
            % c7 = c1*c5    = E*I*G*As*L^2
            % c8 = c1*G*As  = E*I*G*As
            % c9 = rho*A
            % c10= rho*A*Ortsfaktor
            % c11= rho*I
            % c12= c9*L/6   = rho*A*L/6
            % c13= E*A/L
            % c14= G*It/L
            % c15= rho*It*L/6
            
            L = this.Length;
            c = this.c;

            GAs = c(3)/L;
            GAsL = c(3);
            GAsLL = c(5);
            EIy = c(1);
            EIz = c(1);

            N_u =  [-1/L 0 0 0 0 0 1/L 0 0 0 0 0];
            N_phi = [0 0 0 -1/L 0 0 0 0 0 1/L 0 0];

            nenner_y = 1/(L * (12*EIy + GAsLL));
            nenner_z = 1/(L * (12*EIz + GAsLL));
            nenner2_y = 1/(12*EIy + GAsLL);
            nenner2_z = 1/(12*EIz + GAsLL);

            V_v1 = nenner_y * (6*GAs * s^2 - 6*GAsL * s - 12*EIy);
            V_theta1 = nenner_y * (3*GAsL * s^2 - 4*s*(3*EIy + GAsLL) + L*(6*EIy + GAsLL));
            V_v2 = -nenner_y * (6*GAs * s^2 - 6*GAsL * s - 12*EIy);
            V_theta2 = nenner_y * (3*GAsL * s^2 + 2*s*(6*EIy - GAsLL) - 6*EIy*L);

            N_v = [0 V_v1 0 0 0 V_theta1 0 V_v2 0 0 0 V_theta2];


            W_w1 = nenner_z * (6*GAs * s^2 - 6*GAsL * s - 12*EIz);
            W_psi1 = -nenner_z * (3*GAsL * s^2 - 4*s*(3*EIz+GAsLL) + L*(6*EIz + GAsLL));
            W_w2 = -nenner_z * (6*GAs * s^2 - 6*GAsL * s - 12*EIz);
            W_psi2 = -nenner_z * (3*GAsL * s^2 + 2*s*(6*EIz - GAsLL) - 6*EIz*L);

            N_w = [0 0 W_w1 0 W_psi1 0 0 0 W_w2 0 W_psi2 0];

            Psi_w1   = nenner_z *6*GAs * (L - 2*s);
            Psi_psi1 = nenner2_z * 6*GAs * s - nenner_z * 4*(3*EIz + GAsLL);
            Psi_w2   = nenner_z * 6*GAs * (2*s - L);
            Psi_psi2 = nenner2_z * 6*GAs * s + nenner_z * 2*(6*EIz - GAsLL);

            N_psi = [0 0 Psi_w1 0 Psi_psi1 0 0 0 Psi_w2 0 Psi_psi2 0];

            Theta_v1 = nenner_y * 6*GAs * (2*s - L);
            Theta_theta1 = nenner2_y * 6*GAs * s - nenner_y * 4*(3*EIy + GAsLL);
            Theta_v2 = nenner_y * 6*GAs * (L - 2*s);
            Theta_theta2 = nenner2_y * 6*GAs * s + nenner_y * 2*(6*EIy - GAsLL);

            N_theta = [0 Theta_v1 0 0 0 Theta_theta1 0 Theta_v2 0 0 0 Theta_theta2];

            N = [N_u; N_v; N_w; N_phi; N_psi; N_theta];
        end
        
        function plot(this, p, u1, u2, col1, col2, plot_options)
            % Zeichnet Timoshenko Balken mit analytischen Ansatzfunktionen mit Stoffparametern c
            % Die Ansatzfunktionen werden dabei an N Zwischenstellen ausgewertet, d.h. N=0 ist zulässig
            % Anfangs- und Endpunkte (x, y, z), sowie Anfangs- und Endverschiebung (lokal!) (u, v, w, phi, psi, theta) sind gegeben.
            % Colormap-Indices für Anfangs-/Endpunkt col1/col2 sind gegeben, Farbe wird linear interpoliert
            
            % centerline = 2;     % Art der Mittellinie (0: keine, 1: dünn, 2: dick mit Endmarkierung)
            % crosssection = 0;   % Querschnitte (0: keine, 1: nach Theorie 1. Ordn, 2: Theorie 2. Ordn, sonst: exakt rotiert)

            % Create the former function arguments
            split = this.split;
            T = this.T;
            c = this.c;
            p1 = p(this.PointsIdx(1),:);
            p2 = p(this.PointsIdx(2),:);
            %L = norm(p1-p2);
            L = this.Length;

            x = (0:L/(split+1):L)';
            p = x/L*p2 + (L-x)/L*p1;
            col = round(x/L*col2 + (L-x)/L*col1);

            % Ausgangslage
            % plot3(p(:,1), p(:,2), p(:,3),'k:');

            u_nodes = [u1;u2];
            u_lok = zeros(6,split+2);
            u_lok(:,1) = u1;
            u_lok(:,split+2) = u2;

            for i = 2:split+1
                N_U = models.beam.StraightBeam.beam_shape_functions(x(i), L, c);
                u_lok(:,i) = N_U * u_nodes;
            end

            u_glob = [T * u_lok(1:3,:); T * u_lok(4:6,:)];

            cmap = colormap;

            %% Mittellinie plotten
            % Einfarbig, dünn
            if (plot_options.centerline)
            %     plot3(p(:,1) + u_glob(1,:)', p(:,2) + u_glob(2,:)', p(:,3) + u_glob(3,:)', 'Color', 0.5 * (cmap(col(1),:)+cmap(col(split+2),:)) )
            % elseif (centerline == 2)
                % Mehrfarbig, dick
                for i=1:split+1
                    plot3(p([i i+1],1) + u_glob(1,[i i+1])', ...
                          p([i i+1],2) + u_glob(2,[i i+1])', ...
                          p([i i+1],3) + u_glob(3,[i i+1])', ...
                          'LineWidth', plot_options.centerline, ...
                          'Color', 0.5 * (cmap(col(i),:)+cmap(col(i+1),:)) );
                end
            end

            if (plot_options.endmarker)
                % Elementgrenzenmarkierungen plotten
                plot3(p(1,1) + u_glob(1,1)', p(1,2) + u_glob(2,1)', ...
                      p(1,3) + u_glob(3,1)', '+', ...
                      'LineWidth', plot_options.endmarker, 'Color', cmap(col(1),:) );
                plot3(p(split+2,1) + u_glob(1,split+2)', ...
                      p(split+2,2) + u_glob(2,split+2)', ...
                      p(split+2,3) + u_glob(3,split+2)', '+', ...
                      'LineWidth', plot_options.endmarker, 'Color', cmap(col(split+2),:) );
            end

            %% Querschnitte plotten
            if (plot_options.crosssection ~= 0)
                N_angle = 18;
                % R_cross = L/10;
                R_cross = 0.6285;
                angle = 0 : (2*pi/N_angle) : 2*pi;

                xx = zeros(1, N_angle+1);
                yy = R_cross * cos(angle);
                zz = R_cross * sin(angle);

                COR_lok = [xx; yy; zz];

                for i = 1:split+2
                    phi = u_lok(4,i);
                    psi = u_lok(5,i);
                    theta = u_lok(6,i);
                    rot_angle = norm([phi psi theta]);    
                    S_rot = [0 -theta psi; theta 0 -phi; -psi phi 0];

                    % Rotationen 1. Ordnung
                    if (plot_options.crosssection == 1)
                        ROT = eye(3) + S_rot;
                    elseif (plot_options.crosssection == 2)
                        % Rotationen 2. Ordnung
                        ROT = eye(3) + S_rot + 0.5*S_rot*S_rot;
                    else
                        % Exakte Rotationen (R = e^(S_rot))
                        if (rot_angle > 1e-9)
                            ROT = eye(3) + sin(rot_angle)/rot_angle * S_rot + (1-cos(rot_angle))/rot_angle^2 * S_rot*S_rot;
                        else
                            ROT = eye(3);
                        end
                    end

                    COR_u = T * ROT*COR_lok;
                    col_act = cmap(col(i),:);
                    plot3(p(i,1) + u_glob(1,i) + COR_u(1,:), p(i,2) + u_glob(2,i)+ COR_u(2,:), p(i,3) + u_glob(3,i) + COR_u(3,:), 'Color', col_act);
                end
            end
        end
    end
    
    methods(Access=protected)
        function initialize(this)
            % Initializes the straight beam element
            
            % Call superclass
            initialize@models.beam.Beam(this);
            
            m = this.Model;
            pi = this.PointsIdx;

            % Längenberechnung
            dx = (m.Points(pi(2), 1) - m.Points(pi(1), 1));
            dy = (m.Points(pi(2), 2) - m.Points(pi(1), 2)); 
            dz = (m.Points(pi(2), 3) - m.Points(pi(1), 3)); 
            l = sqrt( dx^2 + dy^2 + dz^2 );
            this.Length = l;

            % Lokales Koordinatensystem
            e_x = [dx / l; dy / l; dz / l];
            e_y = [-e_x(2) e_x(1) 0]';
            if norm(e_y) == 0
                e_y = [0 e_x(3) -e_x(2)]';
            end
            e_z = [e_x(2)*e_y(3) - e_x(3)*e_y(2);
                   e_x(3)*e_y(1) - e_x(1)*e_y(3);
                   e_x(1)*e_y(2) - e_x(2)*e_y(1)];
            e_y = e_y / norm(e_y);
            e_z = e_z / norm(e_z);
            this.T = [e_x e_y e_z];

            % Effektive Konstanten
            %   <rho>	<A>     <E>     <Iy/In>     <Iz/Ib>     <It>        <G>     <k>	<c_th>	<kappa>	<alpha>
            %   1       2       3       4           5           6           7       8   9       10      11
            material = this.Material;
            c(1) = material(3) * material(4);                   % c1            = E*I
            c(2) = c(1)^2;                                      % c2 = c1^2     = (E*I)^2
            c(3) = material(8)*material(7)*material(2)*l;       % c3            = G*As*L
            c(4) = c(3)^2;                                      % c4 = c3^2     = (G*As*L)^2
            c(5) = c(3) * l;                                    % c5 = c3*L     = G*As*L^2
            c(6) = c(5)^2;                                      % c6 = c5^2     = (G*As*L^2)^2
            c(7) = c(1) * c(5);                                 % c7 = c1*c5    = E*I*G*As*L^2
            c(8) = c(1) * material(8)*material(7)*material(2);  % c8 = c1*G*As  = E*I*G*As
            c(9) = material(1) * material(2);                   % c9 = rho*A
            c(10) = material(1) * material(2) * 9.81;           %  q = rho*A*Ortsfaktor
            c(11) = material(1) * material(4);                  % c11= rho*I
            c(12) = c(9) * l / 6;                               % c12= c9*L/6   = rho*A*L/6
            c(13) = material(3) * material(2) / l;              % c13= E*A/L
            c(14) = material(7)*material(6) / l;                % c14= G*It/L
            c(15) = material(1)*material(6)*l / 6;              % c15= rho*It*L/6
            this.c = c;
        end
    end
    
    methods(Static)
        function N = beam_shape_functions(s, L, c)
            % Wertet (analytische) Basisfunktion für einen geraden Balken der Länge L und den Stoffkonstanten c an der Stelle s aus.
            % c1            = E*I
            % c2 = c1^2     = (E*I)^2
            % c3            = G*As*L
            % c4 = c3^2     = (G*As*L)^2
            % c5 = c3*L     = G*As*L^2
            % c6 = c5^2     = (G*As*L^2)^2
            % c7 = c1*c5    = E*I*G*As*L^2
            % c8 = c1*G*As  = E*I*G*As
            % c9 = rho*A
            % c10= rho*A*Ortsfaktor
            % c11= rho*I
            % c12= c9*L/6   = rho*A*L/6
            % c13= E*A/L
            % c14= G*It/L
            % c15= rho*It*L/6

            GAs = c(3)/L;
            GAsL = c(3);
            GAsLL = c(5);
            EIy = c(1);
            EIz = c(1);

            xsi = s/L;
            N_u =  [1-xsi 0 0 0 0 0 xsi 0 0 0 0 0];
            N_phi = [0 0 0 1-xsi 0 0 0 0 0 xsi 0 0];

            nenner_y = 1/(L * (12*EIy + GAsLL));
            nenner_z = 1/(L * (12*EIz + GAsLL));
            nenner2_y = 1/(12*EIy + GAsLL);
            nenner2_z = 1/(12*EIz + GAsLL);

            V_v1 = nenner_y * (2*GAs * s^3 - 3*GAsL * s^2 - 12*EIy * s) + 1;
            V_theta1 = nenner_y * s * (GAsL * s^2 - 2*s*(3*EIy + GAsLL) + L*(6*EIy + GAsLL));
            V_v2 = -nenner_y * s * (2*GAs * s^2 - 3*GAsL * s - 12*EIy);
            V_theta2 = nenner_y * s * (GAsL * s^2 + s*(6*EIy - GAsLL) - 6*EIy*L);

            N_v = [0 V_v1 0 0 0 V_theta1 0 V_v2 0 0 0 V_theta2];


            W_w1 = nenner_z * (2*GAs * s^3 - 3*GAsL * s^2 - 12*EIz * s) + 1;
            W_psi1 = -nenner_z * s * (GAsL * s^2 - 2*s*(3*EIz+GAsLL) + L*(6*EIz + GAsLL));
            W_w2 = -nenner_z * s * (2*GAs * s^2 - 3*GAsL * s - 12*EIz);
            W_psi2 = -nenner_z * s *(GAsL * s^2 + s*(6*EIz - GAsLL) - 6*EIz*L);

            N_w = [0 0 W_w1 0 W_psi1 0 0 0 W_w2 0 W_psi2 0];

            Psi_w1   = nenner_z *6*GAs * (L*s - s^2);
            Psi_psi1 = nenner2_z * 3*GAs * s^2 - nenner_z * 4*(3*EIz + GAsLL) * s + 1;
            Psi_w2   = nenner_z * 6*GAs * (s^2 - L*s);
            Psi_psi2 = nenner2_z * 3*GAs * s^2 + nenner_z * 2*(6*EIz - GAsLL) * s;

            N_psi = [0 0 Psi_w1 0 Psi_psi1 0 0 0 Psi_w2 0 Psi_psi2 0];

            Theta_v1 = nenner_y * 6*GAs * (s^2 - L*s);
            Theta_theta1 = nenner2_y * 3*GAs * s^2 - nenner_y * 4*(3*EIy + GAsLL) * s + 1;
            Theta_v2 = nenner_y * 6*GAs * (L*s - s^2);
            Theta_theta2 = nenner2_y * 3*GAs * s^2 + nenner_y * 2*(6*EIy - GAsLL) * s;

            N_theta = [0 Theta_v1 0 0 0 Theta_theta1 0 Theta_v2 0 0 0 Theta_theta2];

            N = [N_u; N_v; N_w; N_phi; N_psi; N_theta];
        end    
    end
end