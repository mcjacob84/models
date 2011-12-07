function plotSingle(this, t, u)
% plot_single: Plots a single beam configuration for given time and field
% data.
%
% Parameters:
% t: The times `t` @type double
% u: The displacement and heat values at `t` @type colvec
%
% @author Daniel Wirtz @date 2011-09-28
%
% @new{0,5,dw,2011-09-28} Added this function.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

% nBeams = numel(this.Beams);
% 
% title_string = sprintf('Balkenwerk zur Zeit t = %2.2f (Verschiebungen um Faktor %2.0f vergrößert)', t, this.PlotFactor);
% title(title_string, 'FontSize', 12, 'FontWeight', 'bold');
% 
% %% Coloring
% col = this.ColorMap;
% if isempty(this.cbh)
%     colormap(col);
%     this.cbh = colorbar;
%     xlabel(this.cbh, 'Temperatur', 'FontSize', 12);
% end
% % Extract temperatures
% temps = u(4:4:end);
% % For real-time plots: Adopt new heat scale if none set yet
% if (max(temps) > this.maxTemp) || (min(temps) < this.minTemp)
%     this.minTemp = min(temps);
%     this.maxTemp = max(temps);
%     caxis([this.minTemp this.maxTemp]);
% end
% col_index = fix((temps - this.minTemp) / (this.maxTemp - this.minTemp + eps) * (size(col,1)-1)) + 1;
% 
% cla;
% hold on;
% for i = 1:nBeams
%     p = this.Beams(i).p(:);
%     % Ausgangslage
%     plot3(this.Points(p,1), this.Points(p,2), this.Points(p,3), 'k:' );
%     
%     % (End)Indizes der Anfangs- und Endpunkte des Elements im globalen
%     % Verschiebungsvektor
%     indices = 4*(this.data.knot_index(p)-1)+1;
%     % (Verstärkte) Verschiebung der der Anfangs- u Endpunkte
%     u_p = this.PlotFactor * [u(indices) u(indices+1) u(indices+2)];
%     plot3(this.Points(p,1) + u_p(:,1), this.Points(p,2) + u_p(:,2),...
%         this.Points(p,3) + u_p(:,3), '-+', 'LineWidth', ...
%         this.BeamLineWidth, 'Color', col(col_index(i),:));
% end        
% hold off;

    % Plot Options
    plot_options.multiplier = 100;                % Vergrößerung der Verschiebungen
    plot_options.figure = 11;                    % Nummer der Figure, in der geplottet wird
    plot_options.colorbar = 'Normalkraft';       % Titel der Colorbar und gleichzeitig zu visualisierende Größe
                                                % Temperatur, Normalkraft, Querkraft y, Querkraft z, Torsionsmoment, 
                                                % Biegemoment y, Biegemoment z,
                                                % Gesamtquerkraft, Gesamtbiegemoment
    plot_options.axis = 'auto';
    % Testszenario
    % plot_options.axis = [-17 17 -17 2 -2 17];   % Grenzen der Achsen des Plots
    % Rohrleitungen y
    % plot_options.axis = [-2 2 -1 11 -4 4];       % Grenzen der Achsen des Plots
    % Rohrleitungen z
    % plot_options.axis = [-2 2 -4 4 -1 11];       % Grenzen der Achsen des Plots
    % Simpel1
    % plot_options.axis = [-2 17 -1 1 -5 5];       % Grenzen der Achsen des Plots
    % "Spazierstock"-Testcase (Simpel2)
    % plot_options.axis = [-17 26 -1 1 -15 10];       % Grenzen der Achsen des Plots

    plot_options.centerline = 3;                % Linienstärke der Mittellinie (0: keine, 1: dünn, 2: dick mit Endmarkierung)
    plot_options.endmarker = 3;                 % Linienstärke der Endmarker
    plot_options.crosssection = 0;              % Querschnitte (0: keine, 1: nach Theorie 1. Ordn, 2: Theorie 2. Ordn, sonst: exakt rotiert)
    plot_options.ref_config = 1;                % Plotten der Ausgangslage (0: nein, 1: ja)

    RO = this.RO;
    KR = this.KR;
    FH = this.FH;
    p = this.Points;
    video = false;
    knoten_index = this.data.knot_index;

    % split_factor_KR = 15;  % Wie oft wird ein Viertelkreis zerlegt
    split_factor_FH = 5;

    num_elem_RO = size(RO, 2);
    num_elem_KR = size(KR, 2);
    num_elem_FH = size(FH, 2);

    h = figure(plot_options.figure);
    clf
    axis equal;
    % grid on
    % axis([0 11 -5.5 5.5]);
    % Portal
    % axis([0 22 0 1 -2 14 ])
    % Rohrleitungen.txt
    % axis([-2 2 -2 2 0 11])
    % Spirale.txt
    % axis([-5 5 -5 5 -5 1])
    % Testszenario.txt
    axis(plot_options.axis)
    set(h, 'Color', 'White')
    if (video)
        set(gcf, 'Renderer', 'zbuffer');
    end
    hold on;

    if (t == -1)
        title_string = sprintf('Statische Lösung (Verschiebungen um Faktor %2.1f vergrößert)', plot_options.multiplier);
    else
        title_string = sprintf('Balkenwerk zur Zeit t = %2.2f (Verschiebungen um Faktor %2.0f vergrößert)', t, plot_options.multiplier);
    end
    title(title_string, 'FontSize', 15, 'FontWeight', 'bold');



    %% Ausgangslage
    if (plot_options.ref_config)
        for i = 1:num_elem_RO
            plot3( p(RO(i).p(:),1), p(RO(i).p(:),2), p(RO(i).p(:),3), 'k:', 'LineWidth', 1 );
        end
        for i = 1:num_elem_KR
        %     split = max(fix(split_factor_KR / (0.5*pi/KR(i).angle)), 1);
            s = 0 : KR(i).angle/(KR(i).split+1) : KR(i).angle;
            % Parametrisierung des Viertelkreises in lokalen Koords
            x = KR(i).R * cos(s);
            y = KR(i).R * sin(s);
            z = 0*s;
            % Umrechnung in globale Koords und Verschiebung um globale Koords des lokalen Ursprungs KR(i).pc
            COR = KR(i).T * [x; y; z];
            plot3( p(KR(i).pc, 1) + COR(1,:), p(KR(i).pc, 2) + COR(2,:), p(KR(i).pc, 3) + COR(3,:), 'k:' );
        end
    %     for i = 1:num_elem_FH
    %         s = 0 : FH(i).l/(4*split_factor_FH) : FH(i).l;
    %         % Parametrisierung des Viertelkreises in lokalen Koords
    %         y = (FH(i).l/(2*split_factor_FH)) * sin( 2*pi*s / (FH(i).l/split_factor_FH) );
    %         x = s;
    %         z = 0*s;
    %         % Umrechnung in globale Koords und Verschiebung um globale Koords des
    %         % lokalen Ursprungs KR(i).pc
    %         COR = FH(i).T * [x; y; z];
    %         plot3( p(FH(i).p(1), 1) + COR(1,:), p(FH(i).p(1), 2) + COR(2,:), p(FH(i).p(1), 3) + COR(3,:), 'k:', 'LineWidth', 1 );
    %     end
    end

    %% Berechnung der zu visualisierenden Größen
    % VARIANTE 1: Visualisierung nach Normalkräften
    % -------------------------------------------------------------------------
    % N = zeros(num_elem_RO + num_elem_KR, 1);
    % offset = 0;
    % for i = 1:num_elem_RO
    % 
    %     T_block = RO(i).T;
    %     
    %     u1 = [u(7*knoten_index(RO(i).p(1))-6) u(7*knoten_index(RO(i).p(1))-5) u(7*knoten_index(RO(i).p(1))-4)]';
    %     % t1 = [u(6*knoten_index(e(i,1))-2) u(6*knoten_index(e(i,1))-1) u(6*knoten_index(e(i,1)))]';
    %     u2 = [u(7*knoten_index(RO(i).p(2))-6) u(7*knoten_index(RO(i).p(2))-5) u(7*knoten_index(RO(i).p(2))-4)]';
    %     % t2 = [u(6*knoten_index(e(i,2))-2) u(6*knoten_index(e(i,2))-1) u(6*knoten_index(e(i,2)))]';
    %     
    %     u1_lok = T_block' * u1;
    %     u2_lok = T_block' * u2;
    %     % t1_lok = T_block' * t1;
    %     % t2_lok = T_block' * t2;
    %     
    %     N(i) = (u1_lok(1) - u2_lok(1)) / RO(i).l;

    % 
    %     offset = offset + 1;
    % end
    % 
    % for i = 1:num_elem_KR
    % 
    %     indices = 7*knoten_index(KR(i).p(:));       
    %     % Verschiebung der der Anfangs- u Endpunkte
    %     u_elem = [u(indices-6) u(indices-5) u(indices-4) u(indices-3) u(indices-2) u(indices-1)];
    %    
    %     T_block_0 = KR(i).T * KR(i).Fren(0);
    %     T_block_L = KR(i).T * KR(i).Fren(KR(i).angle);    
    %         
    %     u0_lok = T_block_0' * u_elem(1,1:3)';
    %     t0_lok = T_block_0' * u_elem(1,4:6)';
    %     uL_lok = T_block_L' * u_elem(2,1:3)';
    %     tL_lok = T_block_L' * u_elem(2,4:6)';
    %     
    %     x = [u0_lok; t0_lok; uL_lok; tL_lok];
    %     
    %     [Ns, B] = circle_shape_functions(KR(i).R, KR(i).R * KR(i).angle, 0);
    %     
    %     Q = B(7:9,:)*x;
    %     
    %     N(offset + i) = Q(1);  
    % 
    % end
    %
    % N_min = min(N);
    % N_max = max(N);


    % VARIANTE 2: Visualisierung nach Temperatur
    % -------------------------------------------------------------------------
    N = zeros(num_elem_RO + num_elem_KR, 1);
    offset = 0;
    for i = 1:num_elem_RO
        % Temperatur der der Anfangs- u Endpunkte
        T1 = u(7*knoten_index(RO(i).p(1)));
        T2 = u(7*knoten_index(RO(i).p(2))); 
        N(i) = 0.5 * (T1 + T2);

        offset = offset + 1;
    end

    for i = 1:num_elem_KR
        indices = 7*knoten_index(KR(i).p(:));       
        % Temperatur der der Anfangs- u Endpunkte
        u_elem = u(indices);

        N(offset + i) = 0.5*(u_elem(1) + u_elem(2));
    end

    % N_min = -.7e6;
    % N_max = .7e6;
    % N_min = -5;
    % N_max = 5;
    % N_min = min(N) - .5
    % N_max = max(N) + .5

    N_min = 0;
    N_max = 2e6;


    %% Colormaps und Farbindexing
    % col = colormap('hot');
    % Eigene Colormap erstellen (num_col Farben)
    % blau -> grün -> rot 
    % num_col = 128;
    % dc = [0:2/num_col:1]';
    % col = [0*dc dc 1-dc; dc 1-dc 0*dc];
    % blau -> grün -> gelb -> rot
    num_col = 128;
    dc = [0:3/num_col:1]';
    col = [0*dc dc 1-dc; dc 0*dc+1 0*dc; 0*dc+1 1-dc 0*dc];
    colormap(col);
    caxis([N_min N_max])

    if (strcmp(plot_options.colorbar,'') == 0)
        handle = colorbar;
        xlabel(handle, plot_options.colorbar, 'FontSize', 12)
    end

    % Den einzelnen Werten von N den richtigen Index in der Colormap zuweisen
    if ( (N_max - N_min) > 1e-5 )
        col_index = fix( (N - N_min) / (N_max - N_min) * (size(col,1)-1)) + 1;
    else
        col_index = 0*N + fix(size(col,1)/2);
    end

    col_index(col_index > num_col) = 128;
    col_index(col_index < 1) = 1;

    val_min = +1e20;
    val_max = -1e20;

    %% Plotten
    offset = 0;
    for i = 1:num_elem_RO
        % (End)Indizes der Anfangs- und Endpunkte des Elements im globalen
        % Verschiebungsvektor
        indices = 7*knoten_index(RO(i).p(:));
        % (Verstärkte) Verschiebung der der Anfangs- u Endpunkte
    %     u_p = plot_options.multiplier * [u(indices-6) u(indices-5) u(indices-4)];

    %     plot3( p(RO(i).p(:),1) + u_p(:,1), p(RO(i).p(:),2) + u_p(:,2), p(RO(i).p(:),3) + u_p(:,3), '-+', 'LineWidth', 3, 'Color', col(col_index(i),:) );
        offset = offset + 1;

        u1_lok = [RO(i).T' * u([indices(1)-6:indices(1)-4]); RO(i).T' * u([indices(1)-3:indices(1)-1])];
        u2_lok = [RO(i).T' * u([indices(2)-6:indices(2)-4]); RO(i).T' * u([indices(2)-3:indices(2)-1])];

         % Ableitungen der Basisfunktionen
        N_prime_1 = this.beam_shape_functions_derivative(0, RO(i).l, RO(i).c);
        N_prime_2 = this.beam_shape_functions_derivative(RO(i).l, RO(i).l, RO(i).c);
        % Ableitungen der Variablen berechnen
        u1_prime_lok = N_prime_1 * [u1_lok; u2_lok];
        u2_prime_lok = N_prime_2 * [u1_lok; u2_lok];   

        if ( strcmpi(plot_options.colorbar, 'Temperatur') )
            % Temperatur
            val1 = u(7*knoten_index(RO(i).p(1)));
            val2 = u(7*knoten_index(RO(i).p(2)));
        elseif ( strcmpi(plot_options.colorbar, 'Normalkraft') )
            % Normalenkraft
            val1 = RO(i).c(13) * RO(i).l * u1_prime_lok(1);
            val2 = RO(i).c(13) * RO(i).l * u2_prime_lok(1);
        elseif ( strcmpi(plot_options.colorbar, 'Querkraft y') )
            % Querkraft y
            val1 = RO(i).c(3) / RO(i).l * (u1_prime_lok(2) - u1_lok(6));
            val2 = RO(i).c(3) / RO(i).l * (u2_prime_lok(2) - u2_lok(6));
        elseif ( strcmpi(plot_options.colorbar, 'Querkraft z') )
            % Querkraft z
            val1 = RO(i).c(3) / RO(i).l * (u1_prime_lok(3) - u1_lok(5));
            val2 = RO(i).c(3) / RO(i).l * (u2_prime_lok(3) - u2_lok(5));
        elseif ( strcmpi(plot_options.colorbar, 'Gesamtquerkraft') )
            % Gesamtquerkraft
            val1 = RO(i).c(3) / RO(i).l * sqrt((u1_prime_lok(2) - u1_lok(6))^2 + (u1_prime_lok(3) - u1_lok(5))^2);
            val2 = RO(i).c(3) / RO(i).l * sqrt((u2_prime_lok(2) - u2_lok(6))^2 + (u2_prime_lok(3) - u2_lok(5))^2);
        elseif ( strcmpi(plot_options.colorbar, 'Torsionsmoment') )
            % Torsionsmoment
            val1 = RO(i).c(14) * RO(i).l * abs(u1_prime_lok(4));
            val2 = RO(i).c(14) * RO(i).l * abs(u2_prime_lok(4));
        elseif ( strcmpi(plot_options.colorbar, 'Biegemoment y') )    
            % Biegemoment y
            val1 = RO(i).c(1) * u1_prime_lok(5);
            val2 = RO(i).c(1) * u2_prime_lok(5);
        elseif ( strcmpi(plot_options.colorbar, 'Biegemoment z') )    
            % Biegemoment z
            val1 = RO(i).c(1) * u1_prime_lok(6);
            val2 = RO(i).c(1) * u2_prime_lok(6);
        elseif ( strcmpi(plot_options.colorbar, 'Gesamtbiegemoment') )
            % Gesamtbiegemoment
            val1 = RO(i).c(1) * sqrt(u1_prime_lok(5)^2 + u1_prime_lok(6)^2);
            val2 = RO(i).c(1) * sqrt(u2_prime_lok(5)^2 + u2_prime_lok(6)^2);
        else
            % Keine Farbe
            val1 = N_min;
            val2 = N_min;
        end

        if (val1 > val_max)
            val_max = val1;
        end
        if (val2 > val_max)
            val_max = val2;
        end
        if (val1 < val_min)
            val_min = val1;
        end
        if (val2 < val_min)
            val_min = val2;
        end

        % Den farbgebenden Werten col1 u col2 aufgrund der Grenzen N_min/max Farbindex zuweisen
        cols = fix( ([val1 val2] - N_min) / (N_max - N_min) * (size(col,1)-1)) + 1;
        cols(cols > num_col) = num_col;
        cols(cols < 1) = 1;
        this.plot_beam(RO(i).split, RO(i).T, RO(i).c, p(RO(i).p(1),:), p(RO(i).p(2),:), plot_options.multiplier*u1_lok, plot_options.multiplier*u2_lok, cols(1), cols(2), plot_options)
    end

    for i = 1:num_elem_KR
        % (End)Indizes der Anfangs- und Endpunkte des Elements im globalen
        % Verschiebungsvektor
        indices = 7*knoten_index(KR(i).p(:));

        % Verschiebung der der Anfangs- u Endpunkte
        u_elem = [u(indices-6) u(indices-5) u(indices-4) u(indices-3) u(indices-2) u(indices-1)];

    %     split = max(fix(split_factor_KR / (0.5*pi/KR(i).angle)), 1);
    %     split = KR(i).split;
    %     s = ( 0 : (KR(i).angle/split) : KR(i).angle );   

    %     x = KR(i).R * cos(s);
    %     y = KR(i).R * sin(s);
    %     z = 0*s;
    %     COR = KR(i).T * [x; y; z];

        u1_lok = KR(i).T_block1' * u_elem(1,1:3)';
        t1_lok = KR(i).T_block1' * u_elem(1,4:6)';
        u2_lok = KR(i).T_block2' * u_elem(2,1:3)';
        t2_lok = KR(i).T_block2' * u_elem(2,4:6)';
        u1 = [u1_lok; t1_lok];
        u2 = [u2_lok; t2_lok];

        % d/ds(u) - theta x e_1 (Querkräfte ohne Stoffkonstanten)
        Q_temp = KR(i).B3*[u1; u2];
        % d/ds(theta) (Momente ohne Stoffkonstanten)
        M_temp = KR(i).B4*[u1; u2];

        if ( strcmpi(plot_options.colorbar, 'Temperatur') )
            % Temperatur
            val1 = u(7*knoten_index(KR(i).p(1)));
            val2 = u(7*knoten_index(KR(i).p(2)));
        % Schnittgrößen über Element konstant!
        elseif ( strcmpi(plot_options.colorbar, 'Normalkraft') )
            % Normalenkraft
            val1 = KR(i).c(3) * Q_temp(1);
            val2 = val1;
        elseif ( strcmpi(plot_options.colorbar, 'Querkraft y') )
            % Querkraft y
            val1 = KR(i).c(4) * Q_temp(2);
            val2 = val1;
        elseif ( strcmpi(plot_options.colorbar, 'Querkraft z') )
            % Querkraft z
            val1 = KR(i).c(4) * Q_temp(3);
            val2 = val1;
        elseif ( strcmpi(plot_options.colorbar, 'Gesamtquerkraft') )
            % Gesamtquerkraft
            val1 = KR(i).c(4) * sqrt(Q_temp(2)^2 + Q_temp(3)^2);
            val2 = val1;
        elseif ( strcmpi(plot_options.colorbar, 'Torsionsmoment') )
            % Torsionsmoment
            val1 = KR(i).c(2) * abs(M_temp(1));
            val2 = val1;
        elseif ( strcmpi(plot_options.colorbar, 'Biegemoment y') )    
            % Biegemoment y
            val1 = KR(i).c(1) * M_temp(2);
            val2 = val1;
        elseif ( strcmpi(plot_options.colorbar, 'Biegemoment z') )    
            % Biegemoment z
            val1 = KR(i).c(1) * M_temp(3);
            val2 = val1;
        elseif ( strcmpi(plot_options.colorbar, 'Gesamtbiegemoment') )
            % Gesamtbiegemoment
            val1 = KR(i).c(1) * sqrt(M_temp(2)^2 + M_temp(3)^2);
            val2 = val1;
        else
            % Keine Farbe
            val1 = N_min;
            val2 = val1;
        end

        if (val1 > val_max)
            val_max = val1;
        end
        if (val2 > val_max)
            val_max = val2;
        end
        if (val1 < val_min)
            val_min = val1;
        end
        if (val2 < val_min)
            val_min = val2;
        end

        cols = fix( ([val1 val2] - N_min) / (N_max - N_min) * (size(col,1)-1)) + 1;
        cols(cols > num_col) = num_col;
        cols(cols < 1) = 1;
        u1 = plot_options.multiplier * u1;
        u2 = plot_options.multiplier * u2;
        this.plot_circle(KR(i).split, KR(i).T, KR(i).T_block1, KR(i).T_block2, KR(i).Fren, KR(i).R, KR(i).angle, KR(i).B, p(KR(i).pc,:), u1, u2, cols(1), cols(2), plot_options)

    %     x = [u0_lok; t0_lok; uL_lok; tL_lok];
    %     x_lok = zeros(6, length(s));
    %     for j = 1:length(s)
    %         Ns = circle_shape_functions(KR(i).R, KR(i).R*s(j), KR(i).B);
    %         x_lok(:,j) = Ns*x;
    %     end
    %     
    %     u_lok = x_lok(1:3,:);
    %     u_p = 0*u_lok;
    %     u_p(:,1) = plot_factor * u_elem(1,1:3)';
    %     u_p(:,length(s)) = plot_factor * u_elem(2,1:3)';
    %     for k = 2:length(s)-1        
    %         u_p(:,k) = plot_factor * KR(i).T * KR(i).Fren(s(k)) * u_lok(:,k);       
    %     end
    %     plot3( p(KR(i).pc, 1) + COR(1,:) + u_p(1,:), p(KR(i).pc, 2) + COR(2,:) + u_p(2,:), p(KR(i).pc, 3) + COR(3,:) + u_p(3,:), 'LineWidth', 3, 'Color', col(col_index(offset + i),:) );
    %     plot3( p(KR(i).p(:),1) + u_p(1,[1 length(s)])', p(KR(i).p(:),2) +  + u_p(2,[1 length(s)])', p(KR(i).p(:),3) + u_p(3,[1 length(s)])', '+', 'LineWidth', 3, 'Color', col(col_index(offset + i),:) );
    % %     plot3( p(KR(i).pc, 1) + COR(1,:) + u_p(1,:), p(KR(i).pc, 2) + COR(2,:) + u_p(2,:), p(KR(i).pc, 3) + COR(3,:) + u_p(3,:), 'LineWidth', 3);

    end

    for i = 1:num_elem_FH

        % Verschiebungsvektor
        indices = 7*knoten_index(FH(i).p(:));

        % Verschiebung der der Anfangs- u Endpunkte
        u_elem = [u(indices-6) u(indices-5) u(indices-4)];

        if (i==5)
            plot3( p(FH(i).p(:), 1) + plot_options.multiplier *u_elem(:,1), p(FH(i).p(:), 2) + plot_options.multiplier *u_elem(:,2), p(FH(i).p(:), 3) + plot_options.multiplier *u_elem(:,3), 'k', 'LineWidth', 3 );
            continue;
        end

    % 	u0_lok = FH(i).T' * u_elem(1,:)';   
    %     uL_lok = FH(i).T' * u_elem(2,:)';
    %     
    %     u_lok = uL_lok * s/FH(i).l + u0_lok * (1 - s/FH(i).l);   

        s = 0 : FH(i).l/(4*split_factor_FH) : FH(i).l;

        u_p = plot_options.multiplier * (u_elem(2,:)' * s/FH(i).l + u_elem(1,:)' * (1 - s/FH(i).l));
        % Parametrisierung des Viertelkreises in lokalen Koords
        y = (FH(i).l/(2*split_factor_FH)) * sin( 2*pi*s / (FH(i).l/split_factor_FH) );
        x = s;
        z = 0*s;
        % Umrechnung in globale Koords und Verschiebung um globale Koords des
        % lokalen Ursprungs KR(i).pc
        COR = FH(i).T * [x; y; z];
        plot3( p(FH(i).p(1), 1) + COR(1,:) + u_p(1,:), p(FH(i).p(1), 2) + COR(2,:) + u_p(2,:), p(FH(i).p(1), 3) + COR(3,:) + u_p(3,:), 'k' );
    end

    % Geschwindigkeitsvektoren plotten
    % for i = 1:length(knoten_index)
    %     if (knoten_index(i))
    %         quiver3(p(i,1)+plot_factor*u(7*knoten_index(i)-6), p(i,2)+plot_factor*u(7*knoten_index(i)-5), p(i,3)+plot_factor*u(7*knoten_index(i)-4), v(7*knoten_index(i)-6), v(7*knoten_index(i)-5), v(7*knoten_index(i)-4), 'Color', [0 0 1])
    %     end
    % end

    % RotMat = viewmtx(37.5, 30);
    % view(RotMat)

    % Spirale, Rohrleitungen
    view(3)

    % view(-37.5-video,30)
    % camva(10)
    % Portal
    % view(0,0)

    % camup([-1 0 1])
    % campos([0 6 0])
    % camtarget([0 0 0])
    hold off;

    val_min
    val_max

    % Frame für Video speichern und anhängen
    if (video)
        frame = getframe(gcf);
    else
        frame = 0;
    end
end