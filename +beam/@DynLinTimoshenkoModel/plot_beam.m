% Zeichnet Timoshenko Balken mit analytischen Ansatzfunktionen mit Stoffparametern c
% Die Ansatzfunktionen werden dabei an N Zwischenstellen ausgewertet, d.h. N=0 ist zulässig
% Anfangs- und Endpunkte (x, y, z), sowie Anfangs- und Endverschiebung (lokal!) (u, v, w, phi, psi, theta) sind gegeben.
% Colormap-Indices für Anfangs-/Endpunkt col1/col2 sind gegeben, Farbe wird linear interpoliert
function plot_beam(this, split, T, c, p1, p2, u1, u2, col1, col2, plot_options)

% centerline = 2;     % Art der Mittellinie (0: keine, 1: dünn, 2: dick mit Endmarkierung)
% crosssection = 0;   % Querschnitte (0: keine, 1: nach Theorie 1. Ordn, 2: Theorie 2. Ordn, sonst: exakt rotiert)

L = norm(p1-p2);

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
        plot3(p([i i+1],1) + u_glob(1,[i i+1])', p([i i+1],2) + u_glob(2,[i i+1])', p([i i+1],3) + u_glob(3,[i i+1])', 'LineWidth', plot_options.centerline, 'Color', 0.5 * (cmap(col(i),:)+cmap(col(i+1),:)) )
    end
end

if (plot_options.endmarker)
    % Elementgrenzenmarkierungen plotten
    plot3(p(1,1) + u_glob(1,1)', p(1,2) + u_glob(2,1)', p(1,3) + u_glob(3,1)', '+', 'LineWidth', plot_options.endmarker, 'Color', cmap(col(1),:) )
    plot3(p(split+2,1) + u_glob(1,split+2)', p(split+2,2) + u_glob(2,split+2)', p(split+2,3) + u_glob(3,split+2)', '+', 'LineWidth', plot_options.endmarker, 'Color', cmap(col(split+2),:) )
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
