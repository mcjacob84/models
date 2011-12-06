% Preprocessing der Rohrleitungsdaten
function preprocess_data(this)
%
% @todo RO_raw mit helper-variable speichern -> sonst wirkt sukzessiver
% aufruf multiplikativ

gravity = [0; 0; -9.81];

% Beams
nBeams = size(this.RO_raw, 1);
this.Beams = struct('p',{});
for i = 1:nBeams
    this.Beams(i).p = [this.RO_raw(i,2) this.RO_raw(i,3)];
end

% Außendurchmesser (m)
ROHR_d_a = 457e-3;
% Wandstärke (m)
ROHR_s = 40e-3;
% Isolierungsdicke (m)
ROHR_iso = 400e-3;
% Manteldicke (m)
ROHR_mantel = 1e-3;
% Dichte des Stahls (kg/m³)(eingelesen!)
ROHR_rho = this.mat(1, 1);
% Dichte der Isolierung (kg/m³)
ROHR_rho_iso = 100;
% Dichte des Mantels (kg/m³)
ROHR_rho_mantel = 7850;
% Dichte des Mediums (kg/m³)
ROHR_rho_med = 20;
% Querkontraktionszahl
ROHR_ny = 0.3;
% E-Modul (N/m²)(eingelesen!)
ROHR_E = this.mat(1, 3);
% Rohrinnendruck (N/m²)
ROHR_p = 50 * 1e5;

% ROHR_E_factor = (1-ROHR_ny) / ( (1+ROHR_ny)*(1-2*ROHR_ny) )

format long
% Innen-/Außendurchmesser des Rohrs
ROHR_r_a = 0.5 * ROHR_d_a;
ROHR_r_i = ROHR_r_a - ROHR_s;
% Querschnittsfläche für Balken
ROHR_A = pi * ( ROHR_r_a^2 - ROHR_r_i^2 );
% Flächenträgheitsmoment für Balken
ROHR_Iy = 0.25 * pi * ( ROHR_r_a^4 - ROHR_r_i^4 );
% Torsionsträgheitsmoment für Balken
ROHR_It = 2 * ROHR_Iy;
% Schubmodul für Balken
ROHR_G = ROHR_E / ( 2*(1+ROHR_ny) );
% Schubkorrekturfaktor für Balken
m_tmp = ROHR_r_i / ROHR_r_a;
ROHR_k = 6*(1+ROHR_ny)*(1+m_tmp^2)^2 / ( (7+6*ROHR_ny)*(1+m_tmp^2)^2 + (20+12*ROHR_ny)*m_tmp^2);

% Berechnung der durch Medium und Dämmung verursachten zusätzlichen Steckenlast
ROHR_q_plus = pi * ( ROHR_r_i^2 * ROHR_rho_med + ( (ROHR_r_a + ROHR_iso)^2 - ROHR_r_a^2 ) * ROHR_rho_iso + ( (ROHR_r_a + ROHR_iso + ROHR_mantel)^2 - (ROHR_r_a + ROHR_iso)^2) * ROHR_rho_mantel );
% ROHR_q_plus =  107.6937961 + 2.23255711 + 31.02416993;
% ROHR_GAs = ROHR_k * ROHR_G * ROHR_A
% ROHR_EI = ROHR_E * ROHR_Iy
% ROHR_q = (ROHR_rho * ROHR_A + ROHR_q_plus) * gravity

% Veränderung von Material 1 gemäß obiger Eingaben
this.mat(1, 1) = ROHR_rho;
this.mat(1, 2) = ROHR_A;
this.mat(1, 3) = ROHR_E;
this.mat(1, 4) = ROHR_Iy;
this.mat(1, 5) = ROHR_Iy;
this.mat(1, 6) = ROHR_It;
this.mat(1, 7) = ROHR_G;
this.mat(1, 8) = ROHR_k;


KR = [];
RO = [];
FH = [];
used_knots_KR = [];
used_knots_RO = [];
used_knots_FH = [];

% Falls nicht alle Knoten verwendet werden
if (size(this.KR_raw, 1))
    used_knots_KR = unique(this.KR_raw(:,2:3));
end
if (size(this.RO_raw, 1))
    used_knots_RO = unique(this.RO_raw(:,2:3));
end
if (size(this.FH_raw, 1))
    used_knots_FH = unique(this.FH_raw(:,2:3));
end
used_knots = union(used_knots_RO, used_knots_KR);
used_knots = union(used_knots, used_knots_FH);

% Neue Nummerierung der Punkte aufstellen (Position in zB globalen
% Verschiebungsvektor)
knot_index = zeros(1,size(this.Points,1));
for index = 1:length(used_knots)
    knot = used_knots(index);
    knot_index(knot) = index;
end

num_knots = length(used_knots);
num_elem_RO = size(this.RO_raw, 1);
num_elem_KR = size(this.KR_raw, 1);
num_elem_FH = size(this.FH_raw, 1);

this.data = struct('num_knots', num_knots, 'knot_index', knot_index, 'num_elem_RO', num_elem_RO, 'num_elem_KR', num_elem_KR, 'num_elem_FH', num_elem_FH);

if (size(this.Supports, 1) > 0)
    this.data.dirichlet = this.Supports(:, 1);
    this.data.dir_data = this.Supports(:, 2:15);
else
    this.data.dirichlet = [];
    this.data.dir_data = [];
end

if (size(this.Loads,1) > 0)
    this.data.neumann = this.Loads(:,1);
    this.data.neu_data = this.Loads(:, 2:8);
else
    this.data.neumann = [];
    this.data.neu_data = [];
end

%% Verarbeitung der geraden Rohrleitungen RO_raw
% RO = struct([]);
gesamtlaenge = 0;
for i = 1:num_elem_RO
    % Materialsatz
    RO(i).Mat = this.RO_raw(i,1);
    
    % Anfangs- und Endpunkt
    RO(i).p = [this.RO_raw(i,2) this.RO_raw(i,3)];
    
    % Längenberechnung
    dx = (this.Points(RO(i).p(2), 1) - this.Points(RO(i).p(1), 1));
    dy = (this.Points(RO(i).p(2), 2) - this.Points(RO(i).p(1), 2)); 
    dz = (this.Points(RO(i).p(2), 3) - this.Points(RO(i).p(1), 3)); 
    RO(i).l = sqrt( dx^2 + dy^2 + dz^2 );
    
gesamtlaenge = gesamtlaenge + RO(i).l;
    
    % Lokales Koordinatensystem
    e_x = [dx / RO(i).l; dy / RO(i).l; dz / RO(i).l];
    e_y = [-e_x(2) e_x(1) 0]';
    if norm(e_y) == 0
        e_y = [0 e_x(3) -e_x(2)]';
    end
    e_z = [e_x(2)*e_y(3) - e_x(3)*e_y(2);
           e_x(3)*e_y(1) - e_x(1)*e_y(3);
           e_x(1)*e_y(2) - e_x(2)*e_y(1)];
    e_y = e_y / norm(e_y);
    e_z = e_z / norm(e_z);
    RO(i).T = [e_x e_y e_z];
    
    % Effektive Konstanten
%   <rho>	<A>     <E>     <Iy/In>     <Iz/Ib>     <It>        <G>     <k>	<c_th>	<kappa>	<alpha>
%   1       2       3       4           5           6           7       8   9       10      11
    mat = RO(i).Mat;
    RO(i).c(1) = this.mat(mat,3) * this.mat(mat,4);                           % c1            = E*I
    RO(i).c(2) = RO(i).c(1)^2;                                  % c2 = c1^2     = (E*I)^2
    RO(i).c(3) = this.mat(mat,8)*this.mat(mat,7)*this.mat(mat,2)*RO(i).l;            % c3            = G*As*L
    RO(i).c(4) = RO(i).c(3)^2;                                  % c4 = c3^2     = (G*As*L)^2
    RO(i).c(5) = RO(i).c(3) * RO(i).l;                          % c5 = c3*L     = G*As*L^2
    RO(i).c(6) = RO(i).c(5)^2;                                  % c6 = c5^2     = (G*As*L^2)^2
    RO(i).c(7) = RO(i).c(1) * RO(i).c(5);                       % c7 = c1*c5    = E*I*G*As*L^2
    RO(i).c(8) = RO(i).c(1) * this.mat(mat,8)*this.mat(mat,7)*this.mat(mat,2);       % c8 = c1*G*As  = E*I*G*As
    RO(i).c(9) = this.mat(mat,1) * this.mat(mat,2);                           % c9 = rho*A
    RO(i).c(10) = this.mat(mat,1) * this.mat(mat,2) * 9.81;                   %  q = rho*A*Ortsfaktor
    RO(i).c(11) = this.mat(mat,1) * this.mat(mat,4);                          % c11= rho*I
    RO(i).c(12) = RO(i).c(9) * RO(i).l / 6;                     % c12= c9*L/6   = rho*A*L/6
    RO(i).c(13) = this.mat(mat,3) * this.mat(mat,2) / RO(i).l;                % c13= E*A/L
    RO(i).c(14) = this.mat(mat,7)*this.mat(mat,6) / RO(i).l;                  % c14= G*It/L
    RO(i).c(15) = this.mat(mat,1)*this.mat(mat,6)*RO(i).l / 6;                % c15= rho*It*L/6
    
    RO(i).q_lok = (RO(i).c(9) + ROHR_q_plus) * (RO(i).T' * gravity);
    RO(i).c_theta = this.mat(mat,9);
    RO(i).kappa = this.mat(mat,10);
    RO(i).alphaA = this.mat(mat,11) * this.mat(mat,2) * this.mat(mat,3);             % alpha*A*E
end

%	KR	<Mat>	<P1>	<P2>	<PCenter>
for i = 1:num_elem_KR
    % Materialsatz
    KR(i).Mat = this.KR_raw(i,1);
    
    % Anfangs- und Endpunkt
    KR(i).p = [this.KR_raw(i,2) this.KR_raw(i,3)];
    KR(i).pc = this.KR_raw(i, 4);
    
    % Lokales Koordinatensystem in globalem entwickeln
    e_x = ( this.Points(KR(i).p(1),:) - this.Points(KR(i).pc,:) )';
    KR(i).R = norm(e_x);
    e_x = e_x / KR(i).R;
    e_y = (this.Points(KR(i).p(2),:) - this.Points(KR(i).pc,:))' / KR(i).R;
    
    % Öffnungswinkel des Kreissegments
    KR(i).angle = acos(e_x'*e_y);
%     e_x'*e_y
%     fac = (pi/2)/acos(e_x'*e_y)

gesamtlaenge = gesamtlaenge + KR(i).R*KR(i).angle;

    e_y = e_y - (e_x'*e_y) * e_x;
    e_y = e_y / norm(e_y);    

    e_z = [e_x(2)*e_y(3) - e_x(3)*e_y(2);
           e_x(3)*e_y(1) - e_x(1)*e_y(3);
           e_x(1)*e_y(2) - e_x(2)*e_y(1)];
    e_z = e_z / norm(e_z);  
    KR(i).T = [e_x e_y e_z];
    
    % Frenet-Basis (in _lokalen_ Koords) als anonyme Funktion
    KR(i).Fren = @(s) [-sin(s) -cos(s) 0; cos(s) -sin(s) 0; 0 0 1];
    
    KR(i).T_block1 = KR(i).T * KR(i).Fren(0);
    KR(i).T_block2 = KR(i).T * KR(i).Fren(KR(i).angle);
    
    % Für die Auswertung der Ansatzfunktionen und Umrechnung der Knotengrößen in Verzerrungen
    B = this.circle_connect_matrix(KR(i).R, KR(i).R * KR(i).angle); 
    KR(i).B = B;
    KR(i).B3 = B(7:9,:);
    KR(i).B4 = B(10:12,:);
    
    % Berechnung des Flexibilitätsfaktors (verringerte Biegesteifigkeit auf Grund von Ovalisierung)
    dm_tmp = ROHR_r_i + ROHR_r_a;
    h_tmp = 4 * KR(i).R * ROHR_s / dm_tmp^2;
    flex_factor = 1.65 / ( h_tmp*(1 + 6*ROHR_p/ROHR_E*(0.5*dm_tmp/ROHR_s)^2*(KR(i).R/ROHR_s)^(1/3)) );
    if (flex_factor < 1)
        flex_factor = 1;
    end
% flex_factor = 1;
    % Effektive Konstanten
%   <rho>	<A>     <E>     <Iy/In>     <Iz/Ib>     <It>        <G>     <k>
%    1       2       3       4           5           6           7       8
    mat = KR(i).Mat;
    KR(i).c(1) = this.mat(mat,3) * this.mat(mat,4) / flex_factor;             % c1 = E*I
    KR(i).c(2) = this.mat(mat,7) * this.mat(mat,6);                           % c2 = G*It
    KR(i).c(3) = this.mat(mat,3) * this.mat(mat,2);                           % c3 = E*A
    KR(i).c(4) = this.mat(mat,8) * this.mat(mat,7) * this.mat(mat,2);                % c4 = G*As
    KR(i).c(5) = this.mat(mat,1) * this.mat(mat,2);                           % c5 = rho*A
    KR(i).c(6) = this.mat(mat,1) * this.mat(mat,4);                           % c6 = rho*I
    KR(i).c(7) = this.mat(mat,1) * this.mat(mat,6);                           % c7 = rho*It
    KR(i).c(8) = this.mat(mat,1);                                      % c8 = rho
    
%     KR(i).q_lok = KR(i).c(5) * (KR(i).T' * gravity);
    KR(i).q_lok = (KR(i).c(5) + ROHR_q_plus)* (KR(i).T' * gravity);
    
    KR(i).c_theta = this.mat(mat,9);
    KR(i).kappa = this.mat(mat,10);
    KR(i).alphaA = this.mat(mat,11) * this.mat(mat,2) * this.mat(mat,3);
    
    % Winkel
    % KR(i).n = this.KR_raw(i,5);
    
end

%	FH	<Mat>	<P1>	<P2>
for i = 1:num_elem_FH
    % Materialsatz
    FH(i).Mat = this.FH_raw(i,1);
    
    % Anfangs- und Endpunkt
    FH(i).p = [this.FH_raw(i,2) this.FH_raw(i,3)];
        
   % Längenberechnung
    dp = (this.Points(FH(i).p(2), :) - this.Points(FH(i).p(1), :));
    
    FH(i).l = norm( dp );    
    
    % Lokales Koordinatensystem
    e_x = dp' / FH(i).l;
    e_y = [-e_x(2) e_x(1) 0]';
    if norm(e_y) == 0
        e_y = [0 e_x(3) -e_x(2)]';
    end
    e_z = [e_x(2)*e_y(3) - e_x(3)*e_y(2);
           e_x(3)*e_y(1) - e_x(1)*e_y(3);
           e_x(1)*e_y(2) - e_x(2)*e_y(1)];
    e_y = e_y / norm(e_y);
    e_z = e_z / norm(e_z);
    FH(i).T = [e_x e_y e_z];
    
%     Effektive Konstanten   
%     <rho>	<A>     <E> 
%       1    2       3
    mat = FH(i).Mat;
    FH(i).c(1) = this.mat(mat,3) * this.mat(mat,2) / FH(i).l;                           % c1 = E*A/L      (Federhärte)
    FH(i).c(2) = this.mat(mat,1) * this.mat(mat,2) * FH(i).l / 6;                       % c2 = rho*A*L/6
    
    FH(i).q_lok = this.mat(mat,1) * this.mat(mat,2) * (FH(i).T' * gravity);
    
    FH(i).c_theta = this.mat(mat,9);
    FH(i).kappa = this.mat(mat,10);
    FH(i).alpha = this.mat(mat,11);
    
end


% Anteil der jeweiligen Elementen an der Gesamtsystemlänge speichern
% Speichern, an wievielen Stellen das Element ausgewertet werden soll, wenn auf das Gesamtsystem M auswertungen kommen sollen
M = 75;
for i = 1:num_elem_RO
    RO(i).split = round(M * RO(i).l / gesamtlaenge); %#ok<*AGROW>
end
for i = 1:num_elem_KR
    KR(i).split = round(M * KR(i).angle * KR(i).R / gesamtlaenge);
end

format short

disp(gesamtlaenge);

this.RO = RO;
this.KR = KR;
this.FH = FH;
