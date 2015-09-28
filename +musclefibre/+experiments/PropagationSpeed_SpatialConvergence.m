% This script performs a spatial convergence test for fixed fibre length
% and decreasing dx resolution. The dx values are chosen orders of
% magnitude around the values from T. Heidlauf's Thesis dx of 0.0052.

%% Clear
clear classes;
m = models.musclefibre.Model('N',2);
if m.System.MotoSarcoLinkIndex ~= 1
    error('Must have MotoSarcoLinkIndex=1 for this experiment.')
end

%% Init
parallel = true;
chunksize = 4; % The maximum single trajectory size in GB
fibrelength = .5; % [cm]
T = 50;
DX = logspace(-1,-4,10);
usenoise = true;
dt = .005;

%% General
minV = -20; %[mV]
buffer = 10; % number of buffer sarcomeres at boundaries, padded to the ones over which velo is measured

% Pad by three sarcomeres to each end
distN = fibrelength./DX;
N = distN + 2*buffer;
nn = length(N);

%% Crunch
base = fullfile(KerMor.App.DataDirectory,'musclefibre','propagationspeed');
tag = sprintf('spconv_T%d_dt%g_noise%d',T,dt,usenoise);
datafile = fullfile(base,['data_' tag '.mat']);

Vms = cell(1,nn);
matlabpool open;
parfor p = 1:nn
    m = models.musclefibre.Model('N',N(p),'SarcoVersion',1,'Noise',usenoise);
    m.T = T;
    m.dt = dt;
    sys = m.System;
    pos = [sys.dm+sys.ds+buffer*sys.dsa+1 % fourth sarco = start
      sys.dm+sys.ds+(N(p)-1-buffer)*sys.dsa+1]; % N-3rd sarco = end
    [~,y] = m.simulate;
    fprintf('Finished run %d with N=%d',p,N(p));
    Vms{p} = y(pos,:);
end
matlabpool close;
save(datafile,'distN','N','Vms');

