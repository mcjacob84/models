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
sv = 1;

file = fullfile(fileparts(mfilename('fullpath')),sprintf('sarcoscaling_v%d.mat',sv));
if exist(file,'file') ~= 2
    m = models.motorunit.Shorten('SarcoVersion',sv,'SPM',true,...
        'DynamicIC',true,'OutputScaling',false);

    m.T = 150;
    m.dt = .1;
    m.EnableTrajectoryCaching = false;

    % Directly set parameter set
    p = [linspace(0,1,60); ones(1,60)*9]; % use max mean current, only need single peak

    n = size(p,2);
    maxvals = zeros(1,n);
    pi = ProcessIndicator('Gathering max forces',n);
    for k=1:n
%     parfor k=1:n
        [t, y] = m.simulate(p(:,k),1);%#ok
        [maxvals(k), pos] = max(y(2,:));
        if pos == length(t)
            error('Time interval too short!')
        end
        pi.step;%#ok
    end
    pi.stop;
    save(file,'p','maxvals');
end

load(file);
pm = PlotManager(false,2,2);

x = p(1,:);
y = 1./maxvals;
h = pm.nextPlot('all');
plot(h, x, y);

% Polynomial fit
reltol = 0.01;
abstol = 0.01;
relerr = inf;
abserr = inf;
deg = 12;
while (relerr > reltol || abserr > abstol) && deg < 20
    c = polyfit(x,y,deg);
    ay = polyval(c,x);
    relerr = max(abs((ay-y)./y));
    abserr = max(abs(ay-y));
    deg = deg+1;
end
h = pm.nextPlot('both');
plot(h,x,y);
hold(h,'on');
ax = 0:.01:1;
ay = polyval(c,ax);
plot(h,ax,ay,'r');

% Error
h = pm.nextPlot('err');
ay = polyval(c,x);
plot(h,x,y-ay,'r',x,(y-ay)./y,'b');
relerr = max(abs((ay-y)./y))
abserr = max(abs(ay-y))
save(file,'c','-APPEND');