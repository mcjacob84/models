% Simulation script for generation of detailed trajectories
% in order to compute action potential propagation speeds and amplitudes
% over time for different motorunits (mean current and fibre type)
%
% See the Processor class for data postprocessing.

%% Clear
clear classes;
%% Init
chunksize = 4; % [GB] The maximum single trajectory size
% The number of sarcomeres over which to measure the propagation speed.
%
% This value has been determined by the
% models.musclefibre.experiments.PropagationSpeed_SpeedConvergence* tests.
%
% for dt=0.005 distN=80 is acceptable, and for double size dt=0.01 "only"
% distN=100 is also okay. As this about halves the storage size, we use
% this variant
distN = 100;
dt = .01;
T = 10000;
usenoise = false;

%% General
minV = -20; %[mV]
buffer = 10; % number of buffer sarcomeres at boundaries, padded to the ones over which velo is measured

% Pad by buffer sarcomeres to each end
N = distN + 2*buffer;
m = models.musclefibre.Model('N',N,'SarcoVersion',1,'Noise',usenoise);
if m.Options.JunctionN ~= 1
    error('Must have MotoSarcoLinkIndex=1 for this experiment.')
end
m.T = T;
% This time-step
m.dt = dt;
t = m.Times;

% Compute time intervals for trajectories of limited size
total_len = length(t);
nchunks = ceil(m.System.NumTotalDofs*8*total_len/(1024^3*chunksize));
interval_idx = floor(linspace(1,total_len,nchunks));
intervals = t(interval_idx);

%% Generate params
% Take a regular mesh and filter by the Motoneuron-Input domain
% See m.System.Params(:).Desired for grid resolution
s = sampling.RandomSampler;
s.Samples = 350;
s.Domain = models.motoneuron.ParamDomain;
mus = s.generateSamples(m);

%% Crunch
base = fullfile(KerMor.App.DataDirectory,'musclefibre_propagationspeed');
tag = sprintf('N%d_T%d_dt%g_noise%d',distN,T,m.dt,usenoise);
thefile = [tag '_data.mat'];
datafile = fullfile(base,thefile);

sys = m.System;
pos = [sys.dm+sys.ds+buffer*sys.dsa+1 % buffer'st sarco = start
      sys.dm+sys.ds+(N-1-buffer)*sys.dsa+1]; % N-buffer'st sarco = end
nmu = size(mus,2);
Vms = cell(1,nmu);
ctimes = zeros(1,nmu);
PCPool.open;
%     for p = 1:nmu
parfor p = 1:nmu
    m = models.musclefibre.Model('N',N,'SarcoVersion',1,'Noise',usenoise);
    m.dt = dt;
    Vm = zeros(2,total_len);
    fprintf('Param %d of %d: Computing %d time intervals\n',p,nmu,length(intervals));
    for iidx = 1:length(intervals)-1
        m.t0 = intervals(iidx);
        m.T = intervals(iidx+1);
        [~,y,ctimes(p)] = m.simulate(mus(:,p));
        Vm(:,interval_idx(iidx):interval_idx(iidx+1)) = y(pos,:);%#ok
        m.System.x0 = dscomponents.ConstInitialValue(y(:,end));
    end
    Vms{p} = Vm;
end
PCPool.close;
save(datafile,'mus','distN','N','t','Vms','ctimes');

