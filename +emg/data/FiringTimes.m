% Script to pre-compute several action potential shapes for faster use
% within the emg model. This is optional, as the correct shapes for the
% selected motor-units can be computed directly for the fibre types used
% in each emg model.
numsamples = 200;
m = models.motoneuron.Model;
% 5000 ms should be enough in terms of precomputed times
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
parfor i = 1:numsamples
    [t,y] = m.simulate(mus(:,i), 1);%#ok
    FT = (y(2,2:end) > 40).*(y(2,1:end-1) <= 40);
    Times{i} = t(logical(FT));
    pi.step;%#ok
end
pi.stop;
save ShapeData Times mus;