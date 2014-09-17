% Script to determine the maximum force created by the sarcomere at maximum
% mean current for each fibre type
%
%
file = fullfile(fileparts(mfilename('fullpath')),'maxforces');

% Dynamic IC: yes, SinglePeak: no, OutScaling: no
m = models.motorunit.Shorten(true, false, false);

m.T = 2000;
m.dt = .1;
m.EnableTrajectoryCaching = false;

% s = m.Sampler;
% s.Samples = 1000;
% s.Seed = 1;
% m.off1_createParamSamples;
% p = m.Data.ParamSamples;
% % Sort samples by fibre type
% [~, idx] = sort(p(1,:));
% p = p(:,idx);
% n = size(p,2);

% Directly set parameter set
n = 100;
p = [linspace(0,1,n); 10*ones(1,n)];

maxvals = zeros(1,n);
maxpos = maxvals;
pi = ProcessIndicator('Gathering max forces',n);
parfor k=1:n
    [t, y] = m.simulate(p(:,k),1);%#ok
    [maxvals(k), pos] = max(y(2,:));
    maxpos(k) = pos;
    if pos == m.T+1
        fprintf('Max at final time for param %d\n',k);
    end
    pi.step;%#ok
end
pi.stop;
save(file,'m','p','maxvals','maxidx','n');

% load(file);
pm = PlotManager;
pm.LeaveOpen = true;

% Load mean current limiting polynomial
s = load(models.motoneuron.Model.FILE_UPPERLIMITPOLY);
eff_mean_current = polyval(s.upperlimit_poly,p(1,:));

h = pm.nextPlot('all','fibre type','values');
plot(h,p(1,:),maxvals,'r',p(1,:),eff_mean_current,'b');

% Polynomial fit
hold(h,'on');
coeff = polyfit(p(1,:),maxvals,3);
ax = 0:.01:1;
ay = polyval(coeff,ax);
plot(h,ax,ay,'r--');
legend('Max activation','Effective mean current','Approx max activation');
pm.done;
save(file,'coeff','-APPEND');