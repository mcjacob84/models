% Script to pre-compute several action potential shapes for faster use
% within the emg model. This is optional, as the correct shapes for the
% selected motor-units can be computed directly for the fibre types used
% in each emg model.
numsamples = 200;
% Use samples with always full activation to have immediate peaks.
mus = [linspace(0,1,numsamples); ones(1,numsamples)*9];
PCPool.open;
for sv = 1:2
    m = models.motorunit.Shorten('SPM',true,'SarcoVersion',sv);
    % Wont need longer as the first peak is usually at about 6ms
    m.T = 100;
    % Fine time resolution to have sufficient detail for fine emg simulations.
    m.dt = .001;
    Shapes = cell(1,numsamples);
    Times = cell(1,numsamples);
    parfor i = 1:numsamples
        [Shapes{i}, Times{i}] = m.getActionPotentialShape(mus(:,i));%#ok
    end
    save(sprintf('ShapeData_v%d.mat',sv),'Shapes','Times','mus');
end
PCPool.close;