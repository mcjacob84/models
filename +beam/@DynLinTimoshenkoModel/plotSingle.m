function plotSingle(model, t, u)
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

nBeams = numel(model.Beams);

title_string = sprintf('Balkenwerk zur Zeit t = %2.2f (Verschiebungen um Faktor %2.0f vergrößert)', t, model.PlotFactor);
title(title_string, 'FontSize', 12, 'FontWeight', 'bold');

%% Coloring
col = model.ColorMap;
if isempty(model.cbh)
    colormap(col);
    model.cbh = colorbar;
    xlabel(model.cbh, 'Temperatur', 'FontSize', 12);
end
% Extract temperatures
temps = u(4:4:end);
% For real-time plots: Adopt new heat scale if none set yet
if (max(temps) > model.maxTemp) || (min(temps) < model.minTemp)
    model.minTemp = min(temps);
    model.maxTemp = max(temps);
    caxis([model.minTemp model.maxTemp]);
end
col_index = fix((temps - model.minTemp) / (model.maxTemp - model.minTemp + eps) * (size(col,1)-1)) + 1;

cla;
hold on;
for i = 1:nBeams
    p = model.Beams(i).p(:);
    % Ausgangslage
    plot3(model.Points(p,1), model.Points(p,2), model.Points(p,3), 'k:' );
    
    % (End)Indizes der Anfangs- und Endpunkte des Elements im globalen
    % Verschiebungsvektor
    indices = 4*(model.System.KnotenIndex(p)-1)+1;
    % (Verstärkte) Verschiebung der der Anfangs- u Endpunkte
    u_p = model.PlotFactor * [u(indices) u(indices+1) u(indices+2)];
    plot3(model.Points(p,1) + u_p(:,1), model.Points(p,2) + u_p(:,2),...
        model.Points(p,3) + u_p(:,3), '-+', 'LineWidth', ...
        model.BeamLineWidth, 'Color', col(col_index(i),:));
end        
hold off;
end