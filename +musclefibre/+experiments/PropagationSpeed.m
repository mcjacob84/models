%% Clear
clear classes;
%% Init
distN = 40; % The number of sarcomeres over which to measure the propagation speed.
T = 500;
minV = -10; %[mV]

% Pad by three sarcomeres to each end
N = distN + 6;

m = models.musclefibre.Model('N',N);
if m.System.MotoSarcoLinkIndex ~= 1
    error('Must have MotoSarcoLinkIndex=1 for this experiment.')
end
m.T = T;
m.dt = m.System.dx/10;

%% Generate params
% Take a regular mesh and filter by the Motoneuron-Input domain
% See m.System.Params(:).Desired for grid resolution
m.Sampler = sampling.GridSampler;
m.Sampler.Domain = models.motoneuron.ParamDomain;
mus = m.Sampler.generateSamples(m);

base = fullfile(KerMor.App.DataDirectory,'musclefibre','propagationspeed');
thefile = sprintf('data_N%d_T%d_dt%g.mat',distN,T,m.dt);
datafile = fullfile(base,thefile);

if exist(datafile,'file') == 2
    load(datafile);
    nmu = size(mus,2);
    sys = m.System;
else
    %% Crunch
    sys = m.System;
    pos = [sys.dm+sys.ds+3*sys.dsa+1 % fourth sarco = start
          sys.dm+sys.ds+(N-1-3)*sys.dsa+1]; % N-3rd sarco = end
    nmu = size(mus,2);
    Vms = cell(1,nmu);
    matlabpool open;
    parfor p = 1:nmu
        [~,y] = m.simulate(mus(:,p)); %#ok<PFBNS>
        % Get V_m values for start and end sarcomeres
        Vms{p} = y(pos,:);
        %plot(t,Vm(1,:),'r',t,Vm(2,:),'b')
    end
    save(datafile,'m','mus','distN','N','minV','Vms');
    matlabpool close;
end

%% Eval
t = m.Times;
V = cell(1,nmu);
P = cell(1,nmu);
pi = ProcessIndicator('Processing',nmu,false);
numi = 100;
tgrid = linspace(0,m.T,numi);
Vinterp = zeros(nmu,numi);
for idx = 1:nmu
    Vm = Vms{idx};
    % Find the peak time indices at junction and end
    p = find(Vm(1,:) > minV);
    peaks_junction = p([1 find(diff(p)>1)+1]);
    p = find(Vm(2,:) > minV);
    peaks_end = p([1 find(diff(p)>1)+1]);
    % Get the time instances the peak was recorded
    peaktimes = [];
    peaktimes(1,:) = t(peaks_junction);
    peaktimes(2,:) = t(peaks_end);
    % Compute velocities
    dist = distN*sys.dx;
    tdiff = diff(peaktimes,1);
    V{idx} = 10*dist./tdiff; % [cm/ms = 0.01m/0.001s = .1m/s]
    P{idx} = peaktimes(1,:);
    Vi = interp1(P{idx},V{idx},tgrid,'cubic');
    %plot(P{idx},V{idx},'r',tgrid,Vi,'b');
    Vinterp(idx,:) = Vi;
    pi.step;
end
pi.stop;

%% Draw
pm = PlotManager;
ax = pm.nextPlot('avgspeed','Average speed (of interpolated velocities) over time for all parameters','time [ms]','speed [m/s]');
plot(ax,tgrid,mean(Vinterp,1));
ax = pm.nextPlot('propspeed','Action potential propagation speed [m/s]','fibre type','mean input current');
tri = delaunay(mus(1,:),mus(2,:));
minV = min(Vinterp(:));
maxV = max(Vinterp(:));
for idx = 1:numi
    if ~ishandle(ax)
        break;
    end
    trisurf(tri,mus(1,:),mus(2,:),Vinterp(:,idx),...
        'FaceColor','interp','EdgeColor','interp','Parent',ax);
    zlim(ax,[minV maxV]);
    title(sprintf('Action potential propagation speed [m/s] at %gs',tgrid(idx)));
    drawnow;
    pause(.1);
end


