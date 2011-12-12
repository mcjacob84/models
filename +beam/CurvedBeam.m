classdef CurvedBeam < models.beam.Beam
    % CurvedBeam:
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
    
    properties
        pc;
        R;
        angle;
        Fren;
        B;
        B3;
        B4;
    end
    
    properties(SetAccess=private)
        T_block1;
        T_block2;
    end
    
    properties(Access=private)
        TM;
    end
    
    methods
        
        function this = CurvedBeam(model, material, pointsidx)
            this = this@models.beam.Beam(model, material, pointsidx(1:2));
            this.pc = pointsidx(3);
            this.initialize;
        end
        
        function M = getLocalMassMatrix(this)
            % Berechnet lokale Steifigkeits- und Massenmatrix eines gekrümmten
            % Timoshenko-Balkens (Viertelkreis) (numerisch mit speziell gewonnenen Ansatzfunktionen)
            %
            % c1 = E*I
            % c2 = G*It
            % c3 = E*A
            % c4 = G*As
            % c5 = rho*A
            % c6 = rho*I
            % c7 = rho*It
                        
            c = this.c;
            
            % Ansatzfunktionen N und Verbindungsmatrix für einen Viertelkreis holen
            
            % Matrizen aufsetzen
            rhoJ = blkdiag(c(7), c(6), c(6));
            
            M = zeros(12);
            
            L = this.Length;
            % Stützstellen für Gaußquadratur
            s = 0.5*L * [ (1-sqrt(3/5)), 1, (1+sqrt(3/5)) ];
            % Gewichte für Gaußquadratur
            w = 0.5*L * 1/9 * [5 8 5];
            
            for i=1:3
                N = this.circle_shape_functions(s(i), this.B);
                N1 = N(1:3,:);
                N2 = N(4:6,:);
                M = M + w(i) * ( c(5)*(N1'*N1) + N2'*rhoJ*N2 );
            end
            
            M = this.TG * M * this.TG';
        end
        
        function K = getLocalStiffnessMatrix(this)
            % Berechnet lokale Steifigkeits- und Massenmatrix eines gekrümmten
            % Timoshenko-Balkens (Viertelkreis) (numerisch mit speziell gewonnenen Ansatzfunktionen)
            %
            % c1 = E*I
            % c2 = G*It
            % c3 = E*A
            % c4 = G*As
            % c5 = rho*A
            % c6 = rho*I
            % c7 = rho*It
                        
            c = this.c;
            % Ansatzfunktionen N und Verbindungsmatrix für einen Viertelkreis holen
            
            % Matrizen aufsetzen
            D = blkdiag(c(3), c(4), c(4));
            E = blkdiag(c(2), c(1), c(1));
            rhoJ = blkdiag(c(7), c(6), c(6));
            
            M = zeros(12);
            
            L = this.Length;
            % Stützstellen für Gaußquadratur
            s = 0.5*L * [ (1-sqrt(3/5)), 1, (1+sqrt(3/5)) ];
            % Gewichte für Gaußquadratur
            w = 0.5*L * 1/9 * [5 8 5];
            
            for i=1:3
                N = this.circle_shape_functions(s(i), this.B);
                N1 = N(1:3,:);
                N2 = N(4:6,:);
                M = M + w(i) * ( c(5)*N1'*N1 + N2'*rhoJ*N2 );
            end
            
            B3 = this.B(7:9,:);
            B4 = this.B(10:12,:);
            
            K = L * (B4' * E * B4 + B3' * D * B3);
            
            K = this.TG * K * this.TG';
        end
        
        function f = getLocalForce(this, gravity)
            % Berechnet lokale Steifigkeits- und Massenmatrix eines gekrümmten
            % Timoshenko-Balkens (Viertelkreis) (numerisch mit speziell gewonnenen Ansatzfunktionen)
            %
            % c1 = E*I
            % c2 = G*It
            % c3 = E*A
            % c4 = G*As
            % c5 = rho*A
            % c6 = rho*I
            % c7 = rho*It
                                    
            f = zeros(12,1);
            
            L = this.Length;
            q_lok = (this.c(5) + this.Material.q_plus)* (this.T' * gravity);
            % Stützstellen für Gaußquadratur
            s = 0.5*L * [ (1-sqrt(3/5)), 1, (1+sqrt(3/5)) ];
            % Gewichte für Gaußquadratur
            w = 0.5*L * 1/9 * [5 8 5];
            for i=1:3
                N = this.circle_shape_functions(s(i), this.B);
                f = f + w(i) * N(1:3,:)'*(this.Fren(s(i)/this.R)'*q_lok);
            end
            
            f = this.TG * f;
        end
        
        function [K, R, U_pot] = getLocalTangentials(this, u)

            % Wertet tangentielle Steifigkeitsmatrix, "Residuumsvektor" (= Vektor der virtuellen internen Energie) und potenzielle Energie eines
            % Timoshenko-Balkens (Kreisbogen mit Öffnungswinkel [angle]) (numerisch mit speziell gewonnenen Ansatzfunktionen) aus

            % c1 = E*I
            % c2 = G*It
            % c3 = E*A
            % c4 = G*As
            % c5 = rho*A
            % c6 = rho*I
            % c7 = rho*It
            
            c = this.c;
            B = this.B;

            % Potenzielle Energie
            U_pot = 0;
            % Residuumsvektor
            R = zeros(12,1);
            % Tangentielle lokale Steifigkeitsmatrix
            K = zeros(12);

            % lokale Verschiebungen des Elements
            u_lok  = this.TG' * u;

            % Stoff-Matrizen aufsetzen
            D = blkdiag(c(3), c(4), c(4));
            E = blkdiag(c(2), c(1), c(1));

            L = this.Length;
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
                N = this.circle_shape_functions(s(i), B);
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
            K = this.TG * K * this.TG';
            R = this.TG * R;
        end
        
        function B = circle_connect_matrix(this)
            C0 = this.circle_shape_functions(0, eye(12));
            % Ehemals: L = this.R * this.angle
            CL = this.circle_shape_functions(this.Length, eye(12));
            B = inv([C0; CL]);
        end
        
        function N = circle_shape_functions(this, s, B)
            % Wertet die Basisfunktionen für ein Kreiselement aus:
            % Auf dem Element werden konstante Strains angenommen => DGLs der Form
            %   x' = Ax + c2
            % müssen gelöst werden. x = (u1, u2, u3, theta1, theta2, theta3). Es ergibt
            % sich eine Lösung der Form:
            %   x(s) = C(s) * c. Wobei C 6x12 und c ein Vektor mit Konstanten der Länge 12
            % ist. (c = [x(0); c2])
            % Um nun Basisfunktionen zu generieren müssen die Lösungen x(s) mit den
            % Werten an den beiden Enden des Elements verbunden werden:
            % v = [x(0); x(L)] = [C(0); C(L)] * c =: inv(B) * c
            
            %       => x(s) = C(s)*B * v =: N(s) * [x(0); x(L)] !
            
            % Da M für jedes Element konstant ist (hängt nur von L ab!), könnte M im
            % Voraus berechnet werden.
            
            % Anmerkung: Die Basisfunktionen sind so gebaut, dass die Ableitung
            % (d/ds u = u' + Gu, etc) eben jene Konstanten sind, die in der zweiten
            % Hälfte von c stehen (c2 s.o.): d/ds N(s)*v = (B*v)(7:12)
            %       (d/ds x(s) = d/ds (C(s))*c = c2 ?)
            
            
            % C = @(s) [
            %     cos(s/R)    sin(s/R)    0   0               0           R-R*cos(s/R)    R*sin(s/R)      R-R*cos(s/R)    0   0                   0                   R*s-R^2*sin(s/R);
            %     -sin(s/R)   cos(s/R)    0   0               0           R*sin(s/R)      R*cos(s/R)-R    R*sin(s/R)      0   0                   0                   R^2-R^2*cos(s/R);
            %     0           0           1   R-R*cos(s/R)    -R*sin(s/R) 0               0               0               s   R*s-R^2*sin(s/R)    R^2*cos(s/R)-R^2    0;
            %     0           0           0   cos(s/R)        sin(s/R)    0               0               0               0   R*sin(s/R)          R-R*cos(s/R)        0;
            %     0           0           0   -sin(s/R)       cos(s/R)    0               0               0               0   R*cos(s/R)-R        R*sin(s/R)          0;
            %     0           0           0   0               0           1               0               0               0   0                   0                   s];
            co = cos(s/this.R);
            si = sin(s/this.R);
            RmRco = this.R-this.R*co;
            R2mR2co = this.R*RmRco;
            Rsi = this.R*si;
            RsmR2si = this.R*s - this.R^2*si;
            Cs =  [
                co	si	0	0       0       RmRco	Rsi     RmRco	0	0           0           RsmR2si;
                -si	co	0   0       0       Rsi     -RmRco	Rsi     0   0           0           R2mR2co;
                0	0	1   RmRco	-Rsi	0       0       0       s   RsmR2si     -R2mR2co    0;
                0	0	0   co      si      0       0       0       0   Rsi         RmRco       0;
                0	0	0   -si     co      0       0       0       0   -RmRco      Rsi         0;
                0	0	0   0       0       1       0       0       0   0           0           s];
            % Cs =  [
            %     co	si	0	0       0       R-R*co	R*si	R-R*co	0	0           0           R*s-R^2*si;
            %     -si	co	0   0       0       R*si	R*co-R	R*si	0   0           0           R^2-R^2*co;
            %     0	0	1   R-R*co	-R*si	0       0       0       s   R*s-R^2*si	R^2*co-R^2	0;
            %     0	0	0   co      si      0       0       0       0   R*si        R-R*co      0;
            %     0	0	0   -si     co      0       0       0       0   R*co-R      R*si        0;
            %     0	0	0   0       0       1       0       0       0   0           0           s];
            N = Cs*B;
        end
        
        function plot(this, p, u1, u2, col1, col2, plot_options)
            % Zeichnet Timoshenko Kreisbogen
            % Die Ansatzfunktionen werden dabei an N Zwischenstellen ausgewertet, d.h. N=0 ist zulässig
            % Anfangs- und Endpunkte (x, y, z), sowie Anfangs- und Endverschiebung (lokal!) (u, v, w, phi, psi, theta) sind gegeben.
            % Colormap-Indices für Anfangs-/Endpunkt col1/col2 sind gegeben, Farbe wird linear interpoliert
            % centerline = 2;     % Art der Mittellinie (0: keine, 1: dünn, 2: dick mit Endmarkierung)
            %
            % crosssection = 0;   % Querschnitte (0: keine, 1: nach Theorie 1. Ordn, 2: Theorie 2. Ordn, sonst: exakt rotiert)

            L = this.Length;
            R = this.R;
            T = this.T;
            N = this.split;
            T1 = this.T_block1;
            T2 = this.T_block2;
            T_Fren = this.Fren;
            B = this.B;
            pc = p(this.pc,:); %#ok<*PROP>

            % Auswertungsstellen
            x = (0:L/(N+1):L)';

            % Ruhelage plotten
            xx = R * cos(x/R);
            yy = R * sin(x/R);
            zz = 0*x;
            COR = T * [xx'; yy'; zz'];
            % plot3( pc(1) + COR(1,:), pc(2) + COR(2,:), pc(3) + COR(3,:), 'k:' );

            u_nodes = [u1;u2];
            u_glob = zeros(6,N+2);
            u_lok = zeros(6,N+2);

            u_glob(:,1) = [T1 * u1(1:3); T1 * u1(4:6)];
            u_glob(:,N+2) = [T2 * u2(1:3); T2 * u2(4:6)];
            u_lok(:,1) = u1;
            u_lok(:,N+2) = u2;

            for i = 2:N+1
                N_U = this.circle_shape_functions(x(i), this.B);
                u_lok(:,i) = N_U * u_nodes;
                T_i = T_Fren(x(i)/R);
                u_glob(:,i) = [T * T_i * u_lok(1:3,i); T * T_i * u_lok(4:6,i)];
            end

            cmap = colormap;
            col = round(x/L*col2 + (L-x)/L*col1);

            %% Mittellinie plotten
            if (plot_options.centerline)
                % Einfarbig, dünn
            %     plot3(pc(1) + COR(1,:) + u_glob(1,:), pc(2) + COR(2,:) + u_glob(2,:), pc(3) + COR(3,:) + u_glob(3,:), 'Color', 0.5 * (cmap(col(1),:)+cmap(col(N+2),:)))
            % elseif (centerline == 2)
                % Mehrfarbig, dick
                for i=1:N+1
                    plot3(pc(1) + COR(1,[i i+1]) + u_glob(1,[i i+1]), pc(2) + COR(2,[i i+1]) + u_glob(2,[i i+1]), pc(3) + COR(3,[i i+1]) + u_glob(3,[i i+1]), 'LineWidth', plot_options.centerline, 'Color', 0.5 * (cmap(col(i),:)+cmap(col(i+1),:)) )
                end
            end

            if (plot_options.endmarker)
                % Elementgrenzenmarkierungen plotten
                plot3( pc(1) + COR(1,1) + u_glob(1,1), pc(2) + COR(2,1) + u_glob(2,1), pc(3) + COR(3,1) + u_glob(3,1), '+', 'LineWidth', plot_options.endmarker, 'Color', cmap(col(1),:) )
                plot3( pc(1) + COR(1,N+2) + u_glob(1,N+2), pc(2) + COR(2,N+2) + u_glob(2,N+2), pc(3) + COR(3,N+2) + u_glob(3,N+2), '+', 'LineWidth', plot_options.endmarker, 'Color', cmap(col(N+2),:) )
            end


            %% Querschnitte plotten
            if (plot_options.crosssection ~= 0)
                N_angle = 18;
                R_cross = 0.6285;
                split_angle = 0 : (2*pi/N_angle) : 2*pi;

                xx = zeros(1, N_angle+1);
                yy = R_cross * cos(split_angle);
                zz = R_cross * sin(split_angle);

                COR_lok = [xx; yy; zz];

                for i = 1:N+2
                    phi = u_lok(4,i);
                    psi = u_lok(5,i);
                    theta = u_lok(6,i);
                    rot_angle = norm([phi psi theta]);    
                    S_rot = [0 -theta psi; theta 0 -phi; -psi phi 0];

                    if (plot_options.crosssection == 1)
                        % Rotationen 1. Ordnung
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

                    COR_u = T * T_Fren(x(i)/R) * ROT*COR_lok;

                    col_act = cmap(col(i),:);
                    plot3(pc(1) + COR(1,i) + u_glob(1,i) + COR_u(1,:), pc(2) + COR(2,i) + u_glob(2,i)+ COR_u(2,:), pc(3) + COR(3,i) + u_glob(3,i) + COR_u(3,:), 'Color', col_act);
                end
            end
        end
    end
    
    methods(Access=protected)
        function initialize(this)
            initialize@models.beam.Beam(this);
            
            m = this.Model;
            % Lokales Koordinatensystem in globalem entwickeln
            e_x = ( m.Points(this.PointsIdx(1),:) - m.Points(this.pc,:) )';
            this.R = norm(e_x);
            e_x = e_x / this.R;
            e_y = (m.Points(this.PointsIdx(2),:) - m.Points(this.pc,:))' / this.R;
            
            % Öffnungswinkel des Kreissegments
            this.angle = acos(e_x'*e_y);
            
            % Länge
            this.Length = this.R*this.angle;
            
            e_y = e_y - (e_x'*e_y) * e_x;
            e_y = e_y / norm(e_y);
            
            e_z = [e_x(2)*e_y(3) - e_x(3)*e_y(2);
                e_x(3)*e_y(1) - e_x(1)*e_y(3);
                e_x(1)*e_y(2) - e_x(2)*e_y(1)];
            e_z = e_z / norm(e_z);
            this.T = [e_x e_y e_z];

            % Frenet-Basis (in _lokalen_ Koords) als anonyme Funktion
            this.Fren = @(s) [-sin(s) -cos(s) 0; cos(s) -sin(s) 0; 0 0 1];
            
            % Globale transformationsmatrix berechnen
            % Transformationsmatrix: natürliche Koords -> glob. Koords
            % Sortierung der Variablen:
            % u0, v0, w0, phi0, psi0, theta0, u1, v1, w1, phi1, psi1, theta1
            %  1   2   3     4     5       6   7   8   9    10    11      12
            this.T_block1 = this.T * this.Fren(0);
            this.T_block2 = this.T * this.Fren(this.angle);
            TG = zeros(12);
            TG([1 2 3], [1 2 3]) = this.T_block1;         % Verschiebung Anfangsknoten
            TG([7 8 9], [7 8 9]) = this.T_block2;         % Verschiebung Endknoten
            TG([4 5 6], [4 5 6]) = this.T_block1;         % Winkel Anfangsknoten
            TG([10 11 12], [10 11 12]) = this.T_block2;   % Winkel Endknoten
            this.TG = TG;
            
            % Für die Auswertung der Ansatzfunktionen und Umrechnung der Knotengrößen in Verzerrungen
            B = this.circle_connect_matrix;
            this.B = B;
            this.B3 = B(7:9,:);
            this.B4 = B(10:12,:);
            
            % Effektive Konstanten
            m = this.Material;
            % Berechnung des Flexibilitätsfaktors (verringerte Biegesteifigkeit auf Grund von Ovalisierung)
            dm_tmp = m.d_a - m.s;
            h_tmp = 4 * this.R * m.s / dm_tmp^2;
            flex_factor = 1.65 / ( h_tmp*(1 + 6*m.p/m.E*(0.5*dm_tmp/m.s)^2*(this.R/m.s)^(1/3)) );
            if (flex_factor < 1)
                flex_factor = 1;
            end
            
            this.c(1) = m.E * m.Iy / flex_factor;   % c1 = E*I
            this.c(2) = m.G * m.It;                 % c2 = G*It
            this.c(3) = m.E * m.A;                  % c3 = E*A
            this.c(4) = m.k * m.G * m.A;            % c4 = G*As
            this.c(5) = m.rho * m.A;                % c5 = rho*A
            this.c(6) = m.rho * m.Iy;               % c6 = rho*I
            this.c(7) = m.rho * m.It;               % c7 = rho*It
            this.c(8) = m.rho;                      % c8 = rho
        end
    end    
end