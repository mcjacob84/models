%% Experiment for determining a good noise input transformation
% The goal is to provide a fibre-type independet scaling of the input mean
% current and noise such that the coefficient of variation of the
% inter-spike intervals is between 15-30%. The test runs assume the
% motoneuron model maps the input mean and noise as follows: 
% 9*((mu(2,:)/9).^mu(3,:))


clear classes;
%%
m = models.motoneuron.Model;
m.TrainingParams = 1:2;
m.T = 2000;

% Manually choose scaling
ns = logspace(log10(.8),log10(2),32);

%%
s = sampling.RandomSampler;
s.Samples = 100;
s.Domain = models.motoneuron.ParamDomain;
mus = s.generateSamples(m);

nns = length(ns);
nmu = length(mus);

%% Compute
covs = zeros(nns,nmu);
freqs = covs;
npeaks = covs;
PCPool.open;
parfor nsidx = 1:nns
    for muidx = 1:nmu
        mu = [mus(:,muidx); ns(nsidx)];%#ok
        [t,y] = m.simulate(mu,1);%#ok
        [peaks, covs(nsidx,muidx), cov_dr, freqs(nsidx,muidx)] = m.analyze(t, y);
        npeaks(nsidx,muidx) = length(peaks);
    end
end
PCPool.close;
save noise_scaling covs freqs mus nmu nns npeaks ns;

%% Or load!
load noise_scaling;

%%
[~,sidx] = sort(mus(1,:));
mus = mus(:,sidx);
covs = covs(:,sidx);
freqs = freqs(:,sidx);
%%
pm = PlotManager(false,4,4);
for k = 1:nns
    cov = covs(k,:);
    npeak = npeaks(k,:);
    ax = pm.nextPlot('',sprintf('ns = %g\nmin/mean(cov)=%g/%g\navg/max peaks=%g/%d',...
        ns(k),min(cov(~isnan(cov))),mean(cov(~isnan(cov))),mean(npeak),max(npeak)),'ft','mc');
    %col = [cov'/max(cov) zeros(nmu,2)];
    
    col = [npeak'/max(npeak) zeros(nmu,2)];
    scatter3(ax,mus(1,:),mus(2,:),cov,15*ones(size(cov)),col,'filled');
    view(ax,[-90 0]);
    %view(ax,[-72 14]);
    zlim(ax,[0 .4]);
end
pm.done;