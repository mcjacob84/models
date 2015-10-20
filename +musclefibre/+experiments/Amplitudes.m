%% Clear
clear classes;%#ok

distN = 100; % The number of sarcomeres over which to measure the propagation speed.
Tfinal = 5000;
dt = .01;
t = 0:dt:Tfinal;
usenoise = false;

%% General
% minV = -20; %[mV]

%% Load data file
base = fullfile(KerMor.App.DataDirectory,'musclefibre','propagationspeed');
tag = sprintf('N%d_T%d_dt%g_noise%d',distN,Tfinal,dt,usenoise);
thefile = [tag '_data.mat'];
load(fullfile(base,thefile));
m = models.musclefibre.Model('N',N,'SarcoVersion',1,'Noise',usenoise);
sys = m.System;
len = distN*sys.dx;

%% Process
p = models.musclefibre.experiments.Processor;
[T,P,invalid] = p.getPeakAmplitudes(t, Vms);
mus(:,invalid) = [];
nmu = size(mus,2);

%% Assemble training data
xi = []; fxi = [];
for k = 1:nmu
    numt = length(T{k});
    xi = [xi [mus(:,k)*ones(1,numt)
              T{k}]];%#ok
    fxi = [fxi P{k}];%#ok
end
ximax = max(xi,[],2);
xi = bsxfun(@times,xi,1./ximax);
atd = data.ApproxTrainData(xi,[],[]);
atd.fxi = fxi;

%% Store results
thefile = [tag '_traindata_amp.mat'];
save(fullfile(base,thefile),'mus','T','P','atd');

%% Compute VKOGA
% Get VKOGA instance
alg = approx.algorithms.VKOGA;
alg.MaxExpansionSize = 600;
alg.UsefPGreedy = false;

nG = 20;
k = [2 3];

% Create configuration for VKOGA
gammas = linspace(.1,2,nG);
comb = Utils.createCombinations(gammas,k);
wc = kernels.config.WendlandConfig('G',comb(1,:),'S',comb(2,:));
wc.Dimension = 3;
ec = kernels.config.ExpansionConfig;
ec.StateConfig = wc;
alg.ExpConfig = ec;

%% Execute computations
kexp = alg.computeApproximation(atd);
%kexp_dbase = kexp.toTranslateBase;

%% Plot the results 1 - visual surface
FunVis2D(kexp, atd);
%% Plot the results 2 - approximation errors
alg.plotSummary;

%% Store results
thefile = [tag '_results_amp.mat'];
save(fullfile(base,thefile),'kexp','ximax','atd','mus');
amp = kexp.toTranslateBase;
save(fullfile(base,'propspeed_amplitudes.mat'),'amp','ximax','-APPEND');

%% Draw
numi = 1000;
pm = PlotManager;
%data = Vinterp;
data = Vpoly;
ax = pm.nextPlot('avgspeed','Average speed (of interpolated velocities) over time for all parameters','time [ms]','speed [m/s]');
plot(ax,tgrid,mean(data,1));

%% Draw 2
ax = pm.nextPlot('propspeed','Action potential propagation speed [m/s]','fibre type','mean input current');
tri = delaunay(mus(1,:),mus(2,:));
minV = min(data(:));
maxV = max(data(:));
zlim(ax,[minV maxV]);
view(ax,[-150 0]);
hold(ax,'on');
trih = [];
for idx = 1:numi
    delete(trih);
    if ~ishandle(ax)
        break;
    end
    trih = trisurf(tri,mus(1,:),mus(2,:),data(:,idx),...
        'FaceColor','interp','EdgeColor','k','Parent',ax);    
    title(sprintf('Action potential propagation speed [m/s] at %gs',tgrid(idx)));
    drawnow;
    pause(.01);        
end
pm.done;

%% Plot discrete signals
pm = PlotManager(false,4,4);
sel1 = find(mus(2,:) > 3);
%sel1 = find(mus(2,:) <= 3);
[~,sel1_sort] = sort(mus(1,sel1),'ascend');
sel1 = sel1(sel1_sort);
nvm = length(sel1);
avgs = zeros(1,nvm);
for idx=1:16%nvm
    k = sel1(idx);
    avgs(idx) = mean(V{k});
    c = polyfit(T{k},V{k},6);
    ax = pm.nextPlot(sprintf('ps_%s_part%d-%d',tag,k,k+15),...
        sprintf('mu=[%g %g], avg v=%g [m/s]',mus(:,k),...
            avgs(idx)),'time','velocity');
    plot(ax,T{k},V{k},'r',tgrid,polyval(c,tgrid),'b');
    %plot(ax,T{k},V{k},'r');
    axis(ax,[0 5000 0 2.1]);
    drawnow;
end

%% Plot plain signals
pm = PlotManager(false,2,3);
pm.ExportDPI = 300;
nvm = length(Vms);
for k=1:nvm
    Vm = Vms{k};
    avgs(k) = mean(V{k});
    ax = pm.nextPlot(sprintf('ps_%s_part%d-%d',tag,k,k+5),...
        sprintf('len=%g, avg v=%g [m/s]',len,...
            avgs(k)),'time','Vm');
    plot(ax,t,Vm(1,:),'r',t,Vm(2,:),'b');
end
%%
pm.savePlots(base,'Format',{'jpg','pdf','fig'},'Close',true);

