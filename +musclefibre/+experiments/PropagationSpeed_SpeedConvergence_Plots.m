% This script performs a convergence test for increasing number of
% sarcomeres at fixed dx resolution.

%% Clear
clear classes;

%% Init
T = 50;
usenoise = true;
dt = .005;

%% Load
base = fullfile(KerMor.App.DataDirectory,'musclefibre','propagationspeed');
tag = sprintf('speedconv_T%d_dt%g_noise%d',T,dt,usenoise);
load(fullfile(base,['data_' tag '.mat']))

%% Process
[T,V,tgrid,Vinterp,Vpoly] = models.musclefibre.experiments.getVelocities(0:dt:T, Vms, measurelengths);

%% Plot
pm = PlotManager(false,2,3);
pm.ExportDPI = 300;
nvm = length(Vms);
avgs = zeros(1,nvm);
for k=1:nvm
    Vm = Vms{k};
    avgs(k) = mean(V{k});
    ax = pm.nextPlot(sprintf('%s_part%d-%d',tag,k,k+5),...
        sprintf('N=%d/len=%g, avg v=%g [m/s]',N(k),measurelengths(k),...
            avgs(k)),'time','Vm');
    plot(ax,t,Vm(1,:),'r',t,Vm(2,:),'b');
end
pm.savePlots(base,'Format',{'jpg','pdf','fig'},'Close',true);
pm = PlotManager;
pm.ExportDPI = 300;
ax = pm.nextPlot(tag,sprintf('SpeedConvergence test %s',...
        tag),'N','average speed [m/s]'); %length [cm]
%plot(ax,measurelengths,avgs);
plot(ax,N-20,avgs);
pm.savePlots(base,'Format',{'jpg','pdf','fig'},'Close',true);