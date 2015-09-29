% Script to pre-compute several action potential shapes for faster use
% within the emg model. This is optional, as the correct shapes for the
% selected motor-units can be computed directly for the fibre types used
% in each emg model.
numsamples = 200;
m = models.motorunit.Shorten('SPM',true,'SarcoVersion',2);
% Wont need longer as the first peak is usually at about 6ms
m.T = 50;
% Fine time resolution to have sufficient detail for fine emg simulations.
m.dt = .001;
% Use samples with always full activation to have immediate peaks.
mus = [linspace(0,1,numsamples); ones(1,numsamples)*9];
Shapes = cell(1,numsamples);
Times = cell(1,numsamples);
pm = PlotManager(false,4,4);
pi = ProcessIndicator('Computing action potential shapes',numsamples,false);
parfor i = 1:numsamples
    [Shapes{i}, Times{i}] = m.getActionPotentialShape(mus(:,i));%#ok
    %ax = pm.nextPlot('dummy',sprintf('Action Potential Shape for mu=[%g,
    %%g]',mus(:,i)),'t [ms]','Vm');
    %plot(ax,Times{i},ShapeData{i});
    pi.step;%#ok
end
pi.stop;
save ShapeData Shapes Times mus