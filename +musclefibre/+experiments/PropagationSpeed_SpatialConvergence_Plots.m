% This script performs a spatial convergence test for fixed fibre length
% and decreasing dx resolution. The dx values are chosen orders of
% magnitude around the values from T. Heidlauf's Thesis dx of 0.0052.

%% Clear
clear classes;

%% Init
T = 50;
DX = logspace(-1,-4,10);
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
pm = PlotManager(false,3,3);
for k=1:length(Vms)
    Vm = Vms{k};
    ax = pm.nextPlot('',sprintf('Length=%g, N=%d/dx=%g',fibrelength,N(k),DX(k)),'time','Vm');
    plot(ax,t,Vm(1,:),'r',t,Vm(2,:),'b');
end
pm.done;

% Gimmick: Computation times plot
% Result from LEAD
% comptimes = [12.7 14.6 19.2 32.7 57.7 307 526 923 2.42e3 6.34e3];
% plot(N,times);
