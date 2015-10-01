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
m = models.musclefibre.Model('N',2);
if m.System.MotoSarcoLinkIndex ~= 1
    error('Must have MotoSarcoLinkIndex=1 for this experiment.')
end

%% Init
fibrelength = .5; % [cm]
T = 50;
DX = logspace(-1,-4,12);
usenoise = true;
dt = .005;

m.T = T;
m.dt = dt;
t = m.Times;

%% General
minV = -20; %[mV]
buffer = 10; % number of buffer sarcomeres at boundaries, padded to the ones over which velo is measured

% Pad by three sarcomeres to each end
distN = round(fibrelength./DX);
N = distN + 2*buffer;
nn = length(N);

%% Crunch
base = fullfile(KerMor.App.DataDirectory,'musclefibre','propagationspeed');
tag = sprintf('spconv_T%d_dt%g_noise%d',T,dt,usenoise);
datafile = fullfile(base,['data_' tag '.mat']);

Vms = cell(1,nn);
if matlabpool('size') == 0
    matlabpool open;
end
parfor p = 1:nn
% for p = 1:nn
    m = models.musclefibre.Model('N',N(p),'dx',DX(p),'SarcoVersion',1,'Noise',usenoise);
    m.T = T;
    m.dt = dt;
    sys = m.System;
    pos = [sys.dm+sys.ds+buffer*sys.dsa+1 % buffer'st sarco = start
      sys.dm+sys.ds+(N(p)-1-buffer)*sys.dsa+1]; % N-buffer'st sarco = end
    [~,y] = m.simulate;
    fprintf('Finished run %d with N=%d',p,N(p));
    Vms{p} = y(pos,:);
end
if matlabpool('size') > 0
    matlabpool close;
end
save(datafile,'N','Vms','t','fibrelength');

% Gimmick: Computation times plot
% Result from LEAD
% comptimes = [12.7 14.6 19.2 32.7 57.7 307 526 923 2.42e3 6.34e3];
% plot(N,times);
