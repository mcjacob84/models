% Zeichnet Timoshenko Kreisbogen
% Die Ansatzfunktionen werden dabei an N Zwischenstellen ausgewertet, d.h. N=0 ist zulässig
% Anfangs- und Endpunkte (x, y, z), sowie Anfangs- und Endverschiebung (lokal!) (u, v, w, phi, psi, theta) sind gegeben.
% Colormap-Indices für Anfangs-/Endpunkt col1/col2 sind gegeben, Farbe wird linear interpoliert
function plot_circle(this, N, T, T1, T2, T_Fren, R, angle, B, pc, u1, u2, col1, col2, plot_options)

% centerline = 2;     % Art der Mittellinie (0: keine, 1: dünn, 2: dick mit Endmarkierung)
% crosssection = 0;   % Querschnitte (0: keine, 1: nach Theorie 1. Ordn, 2: Theorie 2. Ordn, sonst: exakt rotiert)

L = R*angle;

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
    N_U = models.beam.CurvedBeam.circle_shape_functions(R, x(i), B);
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