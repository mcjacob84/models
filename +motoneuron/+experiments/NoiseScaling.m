%%
clear classes;
%%
m = models.motoneuron.Model;
m.TrainingParams = 1:2;
m.T = 300;

% Manually choose scaling
ns = logspace(-1.5,.5,12);

%%
s = sampling.RandomSampler;
s.Samples = 60;
s.Domain = models.motoneuron.ParamDomain;
mus = s.generateSamples(m);

nns = length(ns);
nmu = length(mus);

%%
covs = zeros(nns,nmu);
freqs = zeros(nns,nmu);
PCPool.open;
parfor nsidx = 1:nns
    for muidx = 1:nmu
        mu = [mus(:,muidx); ns(nsidx)];%#ok
        [t,y] = m.simulate(mu,1);%#ok
        [peaks, covs(nsidx,muidx), cov_dr, freqs(nsidx,muidx)] = m.analyze(t, y);
    end
end
PCPool.close;

%%
[~,sidx] = sort(mus(1,:));
mus = mus(:,sidx);
covs = covs(:,sidx);
freqs = freqs(:,sidx);
%%
pm = PlotManager(false,2,2);
for k = 1:nns
    cov = covs(k,:);
    ax = pm.nextPlot('',sprintf('ns = %g, mean(cov)=%g',...
        ns(k),mean(cov(~isnan(cov)))),'ft','mc');
    col = [cov'/max(cov) zeros(nmu,2)];
    scatter3(ax,mus(1,:),mus(2,:),cov,15*ones(size(cov)),col,'filled');
    view(ax,[0 90]);
end
pm.done;