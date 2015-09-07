%% Old version:
% This script runs the motorunit.Model instance without input for 10
% seconds, after which a stable initial condition is found. This is run for
% 50 different mu_1 values (i.e. fibre types), and the resulting values are
% learned by a 7-degree polynomial for each dimension.
% This in turn is used to build a affine-linear initial condition for the
% motorunit.System
%
%% New version:
% In principle the same, but the original shorten IC's are used, a single
% twitch is executed and the value after 2000ms is used as initial
% condition.
%
%% Usage
% The resulting "coeff" matrix is stored in the models.motorunit package
% folder as x0coeff<version>.mat and read upon model construction.

%version = 1;
version = 2;

% Create with no dynamic initial conditions as this script is intended to
% compute them :-)
file = fullfile(fileparts(mfilename('fullpath')),'ic_singletwitch.mat');
filex0 = fullfile(fileparts(mfilename('fullpath')),'..',sprintf('x0coeff%d',version));
if exist(file,'file') == 2
    load(file);
else
    m = models.motorunit.Shorten(version,false,true);
    m.UseNoise = false;
    s = sampling.ManualSampler;
    p = linspace(0,1,200);
    p = [p; ones(size(p))*1];
    s.Samples = p;

    m.off1_createParamSamples;
    m.T = 5000;
    m.dt = .1;

    n = size(p,2);
    maxvals = zeros(1,n);
    maxidx = [];
    endvals = zeros(m.System.f.xDim,n);
    ctimes = zeros(1,n);
    pi = ProcessIndicator('Running %d simulations with one twitch',n/(max(1,matlabpool('size'))),false,n);
    %for k=1:n
    parfor k=1:n
        [t, y, ctimes(k), x] = m.simulate(p(:,k),1);%#ok
        endvals(:,k) = x(:,end);
        pi.step;%#ok
    end
    % Collect all trajectories into the local/main model
    m.Data.TrajectoryData.consolidate(m);
    pi.stop;
    m.save;
    save(file, 'm', 'endvals', 'p', 'ctimes', 'n');
end

%% Fit every dimension by polynomial
deg = 8;
dim = size(endvals,1);
coeff = zeros(dim,deg+1);
err = zeros(dim);
a = dscomponents.AffineInitialValue;
pi = ProcessIndicator('Fitting polynomials',dim);
for k=1:dim
    coeff(k,:) = polyfit(p(1,:),endvals(k,:),deg);
    err(k) = Norm.Linf((endvals(k,:) - polyval(coeff(k,:),p(1,:)))');
    a.addMatrix(sprintf('polyval([%s],mu(1))',sprintf('%g ',coeff(k,:))),full(sparse(k,1,1,dim,1)));
    pi.step;
end
pi.stop;

%% Error checks
endvals_afflin = zeros(size(endvals));
for k=1:n
    endvals_afflin(:,k) = a.evaluate(p(1,k));
end
semilogy(p(1,:),Norm.L2(endvals_afflin-endvals)./Norm.L2(endvals));
title(sprintf('Relative errors of initial values for %d increasing mu_1',n));

save(file, 'a', 'endvals_afflin', 'err', 'coeff', '-APPEND');
save(filex0, 'coeff');