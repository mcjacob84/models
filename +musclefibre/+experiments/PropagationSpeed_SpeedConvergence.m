% This script performs a convergence test for increasing number of
% sarcomeres at fixed dx resolution.
%
%% Results of simulation on LEAD
% Use of 

%% Clear
clear classes;
m = models.musclefibre.Model('N',2);
if m.System.MotoSarcoLinkIndex ~= 1
    error('Must have MotoSarcoLinkIndex=1 for this experiment.')
end

%% Init
T = 50;
usenoise = true;
dt = .01;
m.dt = dt;
m.T = T;
t = m.Times;

%% General
minV = -20; %[mV]
buffer = 10; % number of buffer sarcomeres at boundaries, padded to the ones over which velo is measured

% Pad by three sarcomeres to each end
distN = unique([2:65 round(linspace(2,1000,24))]);
N = distN + 2*buffer;
nn = length(N);
measurelengths = distN*m.System.dx;

%% Crunch
base = fullfile(KerMor.App.DataDirectory,'musclefibre','propagationspeed');
tag = sprintf('speedconv_T%d_dt%g_noise%d',T,dt,usenoise,max(N));
datafile = fullfile(base,['data_' tag '.mat']);

Vms = cell(1,nn);
if matlabpool('size') == 0
    matlabpool open;
end
parfor p = 1:nn
% for p = 1:nn
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
if matlabpool('size') > 0
    matlabpool close;
end
save(datafile,'measurelengths','N','Vms','t');

