function plot(model, t, u)
% plot: 
%
%
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
        
h = figure(1);
%set(h, 'Color', 'White');
axis equal;
xlabel('x'); ylabel('y'); zlabel('z');
%% Plotting bound computations
marginFactor = 1.2;

hlp = reshape(u,4,[]);
mins = min(hlp,[],2);
maxs = max(hlp,[],2);
mins(1:3) = model.PlotFactor*mins(1:3) + min(model.Points,[],1)';
maxs(1:3) = model.PlotFactor*maxs(1:3) + max(model.Points,[],1)';
axis(marginFactor*[mins(1) maxs(1) mins(2) maxs(2) mins(3) maxs(3)]);
% Colors
model.minTemp = mins(4);
model.maxTemp = maxs(4);

% Colorbar stuff
col = model.ColorMap;
colormap(col);
caxis([model.minTemp model.maxTemp]);
model.cbh = colorbar;
xlabel(model.cbh, 'Temperatur', 'FontSize', 12);

for tidx=1:length(t)
    model.plotSingle(t(tidx),u(:,tidx));
    pause(.1);
end