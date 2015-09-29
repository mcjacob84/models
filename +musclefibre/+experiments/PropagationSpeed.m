%% Clear
clear classes;
%% Init
LEAD = false;
if LEAD
    parallel = true; %#ok
    chunksize = 4; % The maximum single trajectory size in GB
else
    parallel = true;
    chunksize = 3; % The maximum single trajectory size in GB
end

% distN = 50; % The number of sarcomeres over which to measure the propagation speed.
% T = 5000;
% usenoise = true;

%% Try with only max/min fibres (see mus = below!)
distN = 300; % The number of sarcomeres over which to measure the propagation speed.
T = 500;
usenoise = false;

% Increase N over fixed dx (longer fibres, speed computation accuracy test)
% Increase N via smaller dx (finer resolution) -> spatial convergence test


%% General
minV = -20; %[mV]
buffer = 10; % number of buffer sarcomeres at boundaries, padded to the ones over which velo is measured

% Pad by three sarcomeres to each end
N = distN + 2*buffer;

m = models.musclefibre.Model('N',N,'SarcoVersion',1,'Noise',usenoise);
if m.System.MotoSarcoLinkIndex ~= 1
    error('Must have MotoSarcoLinkIndex=1 for this experiment.')
end
m.T = T;
m.dt = m.System.dx/10;

% Compute time intervals for trajectories of limited size
nchunks = ceil(m.System.NumTotalDofs*8*length(m.Times)/(1024^3*chunksize));
total_len = length(m.Times);
interval_idx = floor(linspace(1,total_len,nchunks));
intervals = m.Times(interval_idx);

%% Generate params
% Take a regular mesh and filter by the Motoneuron-Input domain
% See m.System.Params(:).Desired for grid resolution
m.Sampler = sampling.GridSampler;
m.Sampler.Domain = models.motoneuron.ParamDomain;
mus = m.Sampler.generateSamples(m);
mus = [0 1; 9 9];

%% Crunch
base = fullfile(KerMor.App.DataDirectory,'musclefibre','propagationspeed');
tag = sprintf('N%d_T%d_dt%g_noise%d',distN,T,m.dt,usenoise);
thefile = ['data_' tag '.mat'];
datafile = fullfile(base,thefile);

if exist(datafile,'file') == 2
    load(datafile);
    nmu = size(mus,2);
    sys = m.System;
else
    sys = m.System;
    pos = [sys.dm+sys.ds+buffer*sys.dsa+1 % fourth sarco = start
          sys.dm+sys.ds+(N-1-buffer)*sys.dsa+1]; % N-3rd sarco = end
    nmu = size(mus,2);
    Vms = cell(1,nmu);
    if parallel
        matlabpool open 2; %#ok
        parfor p = 1:nmu
            m_remote = m;
            orig = m_remote.System.x0;
            Vm = zeros(2,total_len);
            pi = ProcessIndicator('Computing time intervals',length(intervals),false);
            for iidx = 1:length(intervals)-1
                m_remote.t0 = intervals(iidx);
                m_remote.T = intervals(iidx+1);
                [~,y] = m_remote.simulate(mus(:,p));
                Vm(:,interval_idx(iidx):interval_idx(iidx+1)) = y(pos,:);
                m_remote.System.x0 = dscomponents.ConstInitialValue(y(:,end));
                pi.step;
            end
            pi.stop;
            m_remote.System.x0 = orig;            
            Vms{p} = Vm;
        end
        matlabpool close;
    else
        pi = ProcessIndicator('Crunching',nmu*length(intervals),false);
        for p = 1:nmu
            orig = m.System.x0;
            Vm = zeros(2,total_len);
            for iidx = 1:length(intervals)-1
                m.t0 = intervals(iidx);
                m.T = intervals(iidx+1);
                [~,y] = m.simulate(mus(:,p));
                Vm(:,interval_idx(iidx):interval_idx(iidx+1)) = y(pos,:);
                m.System.x0 = dscomponents.ConstInitialValue(y(:,end));
                pi.step;
            end
            m.System.x0 = orig;
            m.t0 = 0;
            m.T = T;
            Vms{p} = Vm;
        end
        pi.stop;
    end
    save(datafile,'m','mus','distN','N','minV','Vms');
end

%% Eval
t = m.Times;
V = cell(1,nmu);
P = cell(1,nmu);
pi = ProcessIndicator('Processing',nmu,false);
numi = 1000;
tgrid = linspace(0,m.T,numi);
Vinterp = zeros(nmu,numi);
Vpoly = zeros(nmu,numi);
maxVs = [];
for idx = 1:nmu
    Vm = Vms{idx};
    
    % Find all locations with negative derivative
    negpos = find(diff(Vm(1,:))<0);
    % Find the positions where the negative derivative starts
    firstpos = [1 find(diff(negpos) > 1)+1];
    % Check that the found locations are peaks and not in the lower noise
    ispeak = Vm(1,negpos(firstpos)) > minV;
    % Remove unwanted
    firstpos(~ispeak) = [];
    peaks_junction = negpos(firstpos);
    % Same for end peaks
    negpos = find(diff(Vm(2,:))<0);
    firstpos = [1 find(diff(negpos) > 1)+1];
    ispeak = Vm(2,negpos(firstpos)) > minV;
    firstpos(~ispeak) = [];
    peaks_end = negpos(firstpos);
    
%     plot(t,Vm(2,:),'r',t(negpos),Vm(2,negpos),'bx')
%     plot(t,Vm(2,:),'r',t(negpos(firstpos)),Vm(2,negpos(firstpos)),'bx')
    % Trivial case
    if length(peaks_junction) > length(peaks_end)
        plot(t,Vm(1,:),'r',t(peaks_junction),Vm(1,peaks_junction),'bx',...
            t,Vm(2,:),'b',t(peaks_end),Vm(2,peaks_end),'rx')
        peaks_junction(end) = [];
        fprintf(2,'Check plot for correct removal of last junction peak!\n');
    end

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
    Vinterp(idx,:) = Vi;
    c = polyfit(P{idx},V{idx},3);
    Vpoly(idx,:) = polyval(c,tgrid);
    pi.step;
end
pi.stop;

%% SVR
% svr = general.regression.ScalarNuSVR;
% svr.nu = .1;
% svr = general.regression.ScalarEpsSVR;
% svr = general.regression.ScalarEpsSVR_SMO;
% svr.Lambda = 0.0005;
% svr.Eps = .01;
% k = kernels.GaussKernel(100);
% %k = kernels.Wendland; k.d = 1; k.k = 3; k.Gamma = 100;
% %kexp = svr.directCompute(k,P{idx},V{idx});
% afx = kexp.evaluate(P{idx});
% tgrid = linspace(0,m.T,1000);
% afx2 = kexp.evaluate(tgrid);
% c = polyfit(P{idx},V{idx},2);
% apol = polyval(c,tgrid);
% plot(P{idx},V{idx},'rx',tgrid,afx2,'b',P{idx},afx,'g',tgrid,apol,'m');

%% Draw
pm = PlotManager;
data = Vinterp;
if size(mus,2) > 2
    ax = pm.nextPlot('avgspeed','Average speed (of interpolated velocities) over time for all parameters','time [ms]','speed [m/s]');
    plot(ax,tgrid,mean(data,1));
    ax = pm.nextPlot('propspeed','Action potential propagation speed [m/s]','fibre type','mean input current');
    tri = delaunay(mus(1,:),mus(2,:));
    minV = min(data(:));
    maxV = max(data(:));
    
    for idx = 1:numi
        if ~ishandle(ax)
            break;
        end
        trisurf(tri,mus(1,:),mus(2,:),data(:,idx),...
            'FaceColor','interp','EdgeColor','interp','Parent',ax);
        zlim(ax,[minV maxV]);
        view(ax,[-150 0]);
        title(sprintf('Action potential propagation speed [m/s] at %gs',tgrid(idx)));
        drawnow;
        pause(.01);        
    end
else
    for k = 1:2
        ax = pm.nextPlot(sprintf('speed_%s_mu%d', tag,k),...
            sprintf('Propagation speed: %s, mu=%s',tag,...
                num2str(mus(:,k))),'t [ms]','velocity [m/s]');
        plot(ax,tgrid,data(k,:),'b',P{k},V{k},'rx');
    end
    pm.savePlots(base,'Format',{'pdf','jpg','fig'});
end
pm.done;

%%
% [~,sidx] = sort(mus(1,:),'ascend');
% [X_,Y_] = meshgrid(mus(1,sidx),tgrid);
% mesh(X_,Y_,Vinterp(sidx,:));
% mesh(Vinterp(sidx,:));
% mesh(Vpoly(sidx,:));


