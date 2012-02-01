function plot(model, t, u)
% plot: Plots the current displacements and heat (u = 4 x numNodes x timesteps matrix)
%
% Parameters:
% t: The times `t_i` @type rowvec
% u: The displacement and heat value matrix, with entry `u_i` in column `i`
% of `u`. @type matrix
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
        
%% Plot Options
% Nummer der Figure, in der geplottet wird
plot_options.figure = 11;                    
plot_options.axis = 'auto';
plot_options.colorbar = 'Normalkraft';

h = figure(plot_options.figure);
axis equal;
axis(plot_options.axis)
% RotMat = viewmtx(37.5, 30);
% view(RotMat)

% Spirale, Rohrleitungen
view(3)

% view(-37.5-video,30)
% camva(10)
% Portal
% view(0,0)
set(h, 'Color', 'White')

if (strcmp(plot_options.colorbar,'') == 0)
    handle = colorbar;
    xlabel(handle, plot_options.colorbar, 'FontSize', 12)
end

ctrl = uicontrol('Tag','btnCancel','Parent',h);
set(ctrl,'String','Cancel','Callback',@(h,e)set(h,'UserData',1));

hold on;

for tidx=1:length(t)
    if get(ctrl,'UserData') == 1
        break;
    end
    model.plotSingle(t(tidx), u(:,tidx), h);
    %pause(.01);
    drawnow;
end