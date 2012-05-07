function [m, r, d] = MathMODMain(dim)
    % Current version works with KerMor 0.6

%     clear classes;
    if nargin == 0
        dim = 1000;
    end
    
    %loc = '/home/dwirtz/aghhome/Material/images';
    loc = 'C:/Users/CreaByte/Documents/Uni/img';
    
    %in = 1;
    
    m = models.mathmod2012.MathMODExperiment(dim);
    m.offlineGenerations;
    r = m.buildReducedModel;
    
    d = tools.EstimatorAnalyzer;
    d.EstimatorIterations = [1 2 5];
    d.EstimatorVersions = [1 1 0 0 1 0 0 1 1];
    d.setModel(r);
    pm = tools.PlotManager(false,2,2);
    pm.FilePrefix = 'mmod_noinput_';
    %d.start([.5; 5; -.2],[],pm);
    d.start([.5; 1; -.2],[],pm);
    pm.done;
%     pm.savePlots(loc,{'fig','pdf','jpg'});
    
%     pm = tools.PlotManager(false,2,2);
%     pm.FilePrefix = 'mmod12';
%     pm.SingleSize = [800 500];
%     % Pure inner dynamics without inputs and initial value -.2
%     tools.ParamSweep(r, [0; 10; -.2], [], 2, 0:1:10, ...
%         pm.nextPlot('mmod_noinput_ps_exppar'));
%     
%     % Pure inner dynamics without inputs and initial value -.2
%     tools.ParamSweep(r, [0; 3; -.2], 1, 1, 0:.1:1, ...
%         pm.nextPlot('mmod_ps_inputpar'));
%     
%     tools.ParamSweep(r, [0; 4; -.2], 2, 1, 0:.1:1, ...
%         pm.nextPlot('mmod_ps_inputpar2'));
%     
%     tools.ParamSweep2D(r, [0; 3; -.2], 1, [1 2], -.5:.2:1.5,-8:1.5:8, ...
%         pm.nextPlot('mmod_ps2d_in_exp'));
%     pm.done;
%     
%     pm.savePlots(loc,{'fig','pdf','jpg'});
end