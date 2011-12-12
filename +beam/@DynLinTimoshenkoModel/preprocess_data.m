% Preprocessing der Rohrleitungsdaten
function preprocess_data(this)
%
% @todo RO_raw mit helper-variable speichern -> sonst wirkt sukzessiver
% aufruf multiplikativ

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

%% Verarbeitung der Randbedingungen (Dirichlet & Neumann)
% Neumann randbedingungen
% @todo mit obiger zuweisung zusammenlegen? wo wirds noch gebraucht?
% (direkt Loads benutzen)
this.f_neum = sparse(7 * num_knots, 1);
for i = 1:length(this.data.neumann)
    if (this.data.knot_index(this.data.neumann(i)) == 0)
        continue;
    end
    index_start = 7*this.data.knot_index(this.data.neumann(i))-6;
    this.f_neum(index_start : index_start+6) = this.f_neum(index_start : index_start+6) + this.data.neu_data(i, :)';
end

% Dirichlet Randbedingungen
nodes = 1:7*num_knots;
% Indices of temperature DoFs
nodes_T = 7 * (1:num_knots);
% Indices of u (space and velocity, 3 each)
nodes_u = setdiff(nodes, nodes_T);

% Aufbau des vollen Dirichlet-Vektors
% Vektor mit Indizes der Freiheitsgrade von u aufbauen
u = zeros(7 * num_knots, 1);
free = 1:7*num_knots;

% Balken
for i = 1:length(this.data.dirichlet)
    % Lokale Liste der Komponenten dieses Knotens, die Dirichlet sind
    dir_knoten_lok = (1:7) .* this.data.dir_data(i,8:14);
    dir_knoten_lok = setdiff(dir_knoten_lok, 0);

    % Wenn der Knoten in der Struktur nicht verwendet wird, auslassen
    if (this.data.knot_index(this.data.dirichlet(i)) == 0)
        continue;
    end

    % Globale Liste der Komponenten dieses Knotens, die Dirichlet sind
    dir_knoten = 7*this.data.knot_index(this.data.dirichlet(i))-6 : 7*this.data.knot_index(this.data.dirichlet(i));
    dir_knoten = dir_knoten .* this.data.dir_data(i,8:14);
    dir_knoten = setdiff(dir_knoten, 0);

    % Dirichlet-Werte in Lösungsvektor eintragen und Liste der Freiheitsgrade verkleinern
    u(dir_knoten) = this.data.dir_data(i, dir_knoten_lok)';
    free = setdiff(free, dir_knoten);
end

%% Vektor der Freiheitsgrade aufteilen in Freiheitsgrade für Balken und Temperatur
% Vektor für Balken (Temperatureinträge entfernen)
free_u = setdiff(free, nodes_T);
if this.withHeat
    % "free" contains all DoFs, both u and T, just extract the T DoFs here
    free_T = setdiff(free, free_u);
else
    % Only consider u DoFs
    free = free_u;
end
this.free = free;

this.dir_u = setdiff(nodes_u, free_u);
if this.withHeat
    this.dir_T = setdiff(nodes_T, free_T);
end
% Store dirichlet values
this.u_dir = u(this.dir_u);

%% Erstellen der Strukturelemente
this.Elements = cell(1,num_elem_RO+num_elem_KR+num_elem_FH);
% Verarbeitung der geraden Rohrleitungen RO_raw
gesamtlaenge = 0;
cnt = 1;
for i = 1:num_elem_RO
    % constructor arguments: model instance, material, pointsindices
    sb = models.beam.StraightBeam(this, this.Materials(this.RO_raw(i,1)), this.RO_raw(i,2:3));
    gesamtlaenge = gesamtlaenge + sb.Length;
    this.Elements{cnt} = sb;
    cnt = cnt+1;
end
% Verarbeitung der krummen Rohrleitungen KR_raw
%	KR	<Mat>	<P1>	<P2>	<PCenter>
for i = 1:num_elem_KR
    cb = models.beam.CurvedBeam(this, this.Materials(this.KR_raw(i,1)), this.KR_raw(i,2:4));
    gesamtlaenge = gesamtlaenge + cb.Length;
    this.Elements{cnt} = cb;
    cnt = cnt+1;
end
% Verarbeitung der Stäbe FH_raw
%	FH	<Mat>	<P1>	<P2>
for i = 1:num_elem_FH
    this.Elements{cnt} = models.beam.Truss(this, this.Materials(this.FH_raw(i,1)), this.FH_raw(i,2:3));
    cnt = cnt+1;
end

%% Anteil der jeweiligen Elementen an der Gesamtsystemlänge speichern
% Speichern, an wievielen Stellen das Element ausgewertet werden soll, wenn auf das Gesamtsystem M auswertungen kommen sollen
% Gilt nur für RO und KR elemente
M = 75;
for i = 1:num_elem_RO+num_elem_KR
    this.Elements{i}.split = round(M * this.Elements{i}.Length / gesamtlaenge);
end

fprintf('Gesamtlänge: %f\n',gesamtlaenge);