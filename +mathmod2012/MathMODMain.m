function [m, r, d] = MathMODMain(dim)
    % Current version works with KerMor 0.6

%     clear classes;
    if nargin == 0
        dim = 240000;
    end
    
    loc = '/home/dwirtz/aghhome/Diss/errorestimation/img';
%     loc = 'C:/Users/CreaByte/Documents/Uni/img';
    
    %in = 1;
    
    m = models.mathmod2012.MathMODExperiment(dim);
    m.offlineGenerations;
    r = m.buildReducedModel;
    
    pm = PlotManager;
    pm.FilePrefix = 'mmod_noinput_';
    pm.NoTitlesOnSave = true;
    pm.AutoTickMarks = false;
    pm.UseFileTypeFolders = false;
    
    %% Estimation plots for fixed param
    d = EstimatorAnalyzer;
    d.EstimatorIterations = [1 2 5];
    d.EstimatorVersions = [1 1 0 0 1 0 0 1 1];
    d.setModel(r);
%     [e,re,ct] = d.compute([.5; 5; -.2], []);
%     d.createPlots(e,re,ct,pm);
%     pm.done;
%     pm.savePlots(loc,'Format','eps');
    
    %% Parameter sweeps
    % Use TD estimator
    r.ErrorEstimator = d.Est(7).Estimator;
    pm.ExportDPI = 300;
    
    % Pure inner dynamics without inputs and initial value -.2
%     ParamSweep(r, [0; 10; -.2], [], 2, 0:1:10, pm);
%     pm.FilePrefix = 'mmod_exppar';
%     pm.savePlots(loc,'Format',{'jpg','fig'},'Close',true);
    
    % Pure inner dynamics without inputs and initial value -.2
    pm.SingleSize = [800 500];
    ParamSweep(r, [0; 3; -.2], 1, 1, 0:.1:1, pm);
    view(40,10);
    pm.FilePrefix = 'mmod_inputpar';
    pm.savePlots(loc,'Format',{'jpg','fig'},'Close',true);
    
    ParamSweep(r, [0; 4; -.2], 2, 1, 0:.1:1, pm);
    pm.FilePrefix = 'mmod_inputpar2';
    view(42,32);
    pm.savePlots(loc,'Format',{'jpg','fig'},'Close',true);
    
    ParamSweep2D(r, [0; 3; -.2], 1, [1 2], -.5:.2:1.5,-8:1.5:8, pm.nextPlot('in_exp_2d'));
    pm.FilePrefix = 'mmod';
    pm.savePlots(loc,'Format',{'jpg','fig'},'Close',true);
end