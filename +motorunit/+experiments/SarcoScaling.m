% Script to create an output force scaling for the Shorten sarcomere model.
% The assumption is to have the same force being generated from each fibre
% type when only one twitch is performed.
%
% This is still questionable!
%
% The resulting polynomial coefficients are stored in the
% models.motorunit.SHSystem.ForceOutputScalingPolyCoeff property and used
% for appropriate scaling of the output signal (depending on the fibre
% type)
file = fullfile(fileparts(mfilename('fullpath')),'sarcoscaling');

% m = models.motorunit.Shorten(true, true);
% 
% m.T = 80;
% m.dt = .1;
% m.EnableTrajectoryCaching = false;
% 
% % s = m.Sampler;
% % s.Samples = 1000;
% % s.Seed = 1;
% % m.off1_createParamSamples;
% % p = m.Data.ParamSamples;
% % % Sort samples by fibre type
% % [~, idx] = sort(p(1,:));
% % p = p(:,idx);
% 
% % Directly set parameter set
% p = [linspace(0,1,60); linspace(2,3,60)];
% 
% n = size(p,2);
% maxvals = zeros(1,n);
% maxidx = [];
% pi = ProcessIndicator('Gathering max forces',n);
% parfor k=1:n
%     [t, y, ct, x] = m.simulate(p(:,k),1);%#ok
%     [maxvals(k), pos] = max(x(59,:));
%     if pos == m.T+1
%         fprintf('Max at final time for param %d\n',k);
%         maxidx = [maxidx k];
%     end
%     pi.step;%#ok
% end
% pi.stop;
% save(file,'m','p','maxvals','maxidx');

load(file);
pm = PlotManager(false,2,2);
tri = delaunay(p(1,:),p(2,:));
h = pm.nextPlot('all');
trisurf(tri, p(1,:),p(2,:),maxvals,'Parent',h, 'FaceColor','interp','EdgeColor','interp');

% m.Sampler.Domain = MotoneuronParamDomain;
[~, idx] = m.Sampler.Domain.filter(p);
% tri = delaunay(p(1,idx),p(2,idx));
% h = pm.nextPlot('valid');
% trisurf(tri, p(1,idx),p(2,idx), maxvals(idx), 'Parent',h, 'FaceColor','interp','EdgeColor','interp');

x = p(1,idx);
y = 1./maxvals(idx);
h = pm.nextPlot('all');
plot(h, x, y);
axis(h,'tight');

% Polynomial fit
h = pm.nextPlot('all');
plot(h,x,y);
hold(h,'on');
c = polyfit(x,y,12);
ax = 0:.01:1;
ay = polyval(c,ax);
plot(h,ax,ay,'r');
axis(h,'tight');
save(file,'c','-APPEND');