% Preprocessing der Rohrleitungsdaten
function preprocess_data(this)
%
% @todo RO_raw mit helper-variable speichern -> sonst wirkt sukzessiver
% aufruf multiplikativ

% Beams
% nBeams = size(this.RO_raw, 1);
% this.Beams = struct('p',{});
% for i = 1:nBeams
%     this.Beams(i).p = [this.RO_raw(i,2) this.RO_raw(i,3)];
% end

% Veränderung von Material 1 gemäß obiger Eingaben
% this.mat(1, 1) = ROHR_rho;
% this.mat(1, 2) = ROHR_A;
% this.mat(1, 3) = ROHR_E;
% this.mat(1, 4) = ROHR_Iy;
% this.mat(1, 5) = ROHR_Iy;
% this.mat(1, 6) = ROHR_It;
% this.mat(1, 7) = ROHR_G;
% this.mat(1, 8) = ROHR_k;

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

this.data = struct('num_knots', num_knots, 'knot_index', knot_index,...
                   'num_elem_RO', num_elem_RO, 'num_elem_KR', num_elem_KR, 'num_elem_FH', num_elem_FH);

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
RO = models.beam.StraightBeam.empty;
gesamtlaenge = 0;
for i = 1:num_elem_RO
    
    material = this.mat(this.RO_raw(i,1),:);
    % constructor arguments: model instance, material, pointsindices
    RO(i) = models.beam.StraightBeam(this, material, this.RO_raw(i,2:3));
        
    gesamtlaenge = gesamtlaenge + RO(i).Length;
end
%% Verarbeitung der krummen Rohrleitungen KR_raw
%	KR	<Mat>	<P1>	<P2>	<PCenter>
KR = models.beam.CurvedBeam.empty;
for i = 1:num_elem_KR
    
    material = this.mat(this.KR_raw(i,1),:);
    KR(i) = models.beam.CurvedBeam(this, material, this.KR_raw(i,2:4));
    
    gesamtlaenge = gesamtlaenge + KR(i).Length;
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
    
    FH(i).q_lok = this.mat(mat,1) * this.mat(mat,2) * (FH(i).T' * this.Gravity);
    
    FH(i).c_theta = this.mat(mat,9);
    FH(i).kappa = this.mat(mat,10);
    FH(i).alpha = this.mat(mat,11);
    
end

% Anteil der jeweiligen Elementen an der Gesamtsystemlänge speichern
% Speichern, an wievielen Stellen das Element ausgewertet werden soll, wenn auf das Gesamtsystem M auswertungen kommen sollen
M = 75;
for i = 1:num_elem_RO
    RO(i).split = round(M * RO(i).Length / gesamtlaenge); %#ok<*AGROW>
end
for i = 1:num_elem_KR
    KR(i).split = round(M * KR(i).angle * KR(i).R / gesamtlaenge);
end

format short

disp(gesamtlaenge);

this.RO = RO;
this.KR = KR;
this.FH = FH;
