% This script performs a spatial convergence test for fixed fibre length
% and decreasing dx resolution. The dx values are chosen orders of
% magnitude around the values from T. Heidlauf's Thesis dx of 0.0052.
%
%% Results from simulation on LEAD
% Use dx ~ 10^-2.5=0.0032 [cm] for sufficient
% spatial resolution. For quick simulations, dx=0.05 is the largest value
% still having a reasonably good resolution.


%% Clear
clear classes;

%% Init
T = 50;
DX = logspace(-1,-4,12);
usenoise = true;
dt = .005;
fibrelength = .5; % [cm]

%% Load
base = fullfile(KerMor.App.DataDirectory,'musclefibre','propagationspeed');
tag = sprintf('spconv_T%d_dt%g_noise%d',T,dt,usenoise);
load(fullfile(base,['data_' tag '.mat']))

%% Process
[Times,V,tgrid,Vinterp,Vpoly] = models.musclefibre.experiments.getVelocities(t, Vms, fibrelength);

%% Plot
pm = PlotManager(false,2,2);
pm.ExportDPI = 300;
nvm = length(Vms);
avgv = zeros(1,nvm);
for k=1:nvm
    Vm = Vms{k};
    ax = pm.nextPlot(sprintf('%s_parts%d-%d',tag,k,k+3),...
        sprintf('Length=%g, N=%d/dx=%g, avg velocity=%g [m/s]',...
        fibrelength,N(k),DX(k),mean(V{k})),'time','Vm');
    plot(ax,t,Vm(1,:),'r',t,Vm(2,:),'b');
    avgv(k) = mean(V{k});
end
pm.savePlots(base,'Format',{'jpg','pdf','fig'},'Close',true);
pm = PlotManager;
pm.ExportDPI = 300;
ax = pm.nextPlot(tag,sprintf('SpatialConvergence test %s',tag),'dx [10^{-x}]','average speed [m/s]');
plot(ax,-log10(DX),avgv);
pm.savePlots(base,'Format',{'jpg','pdf','fig'},'Close',true);

% Gimmick: Computation times plot
% Result from LEAD
% comptimes = [12.7 14.6 19.2 32.7 57.7 307 526 923 2.42e3 6.34e3];
% plot(N,times);
