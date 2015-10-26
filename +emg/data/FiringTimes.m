% Script to pre-compute several firing times for faster use
% within the emg model. This is optional, as the correct firing times for the
% selected motor-units can be computed directly for the fibre types used
% in each emg model.
numsamples = 200;
m = models.motoneuron.Model;
% 10000 ms should be enough in terms of precomputed times
m.T = 10000;
% Fine time resolution to have sufficient detail for fine emg simulations.
m.dt = .1;
% Use samples with always full activation to have immediate peaks.
s = sampling.RandomSampler;
s.Domain = models.motoneuron.ParamDomain;
s.Samples = numsamples;
m.Sampler = s;

m.off1_createParamSamples;
mus = m.Data.ParamSamples;
Times = cell(1,numsamples);
pi = ProcessIndicator('Computing motoneuron firing times',numsamples,false);
PCPool.open;
parfor i = 1:numsamples
    [t,y] = m.simulate(mus(:,i), 1);%#ok
    peaks = m.analyze(t,y);
    Times{i} = t(peaks);
    pi.step;%#ok
end
PCPool.close;
pi.stop;
save FiringTimes Times mus;