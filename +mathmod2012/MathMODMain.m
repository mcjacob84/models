function [m, r, d] = MathMODMain(dim)
    % Current version works with KerMor 0.5

%     clear classes;
    if nargin == 0
        dim = 1000;
    end
    
    loc = '/home/dwirtz/aghhome/Material/images';
    figsize = [800 500];
    
    
    %in = 1;
    
    m = models.mathmod2012.MathMODExperiment(dim);
    m.offlineGenerations;
    r = m.buildReducedModel;
    
%     d = tools.EstimatorAnalyzer;
%     d.EstimatorIterations = [1 2 5];
%     d.EstimatorVersions = [1 1 0 0 1 0 0 1 1];
%     d.SingleFigures = true;
%     d.setModel(r);
%     d.start([.5; 5; -.2],[]);
%     
%     names = {'abserr','relerr','ctimes'};
%     prefix = 'mmod_noinput_';
%     for i = 1:3
%         saveFig(i,fullfile(loc,[prefix names{i}]));
%     end
    
    % Pure inner dynamics without inputs and initial value -.2
    fig = tools.ParamSweep(r, [0; 10; -.2], [], 2, 0:1:10);
    saveFig(fig, fullfile(loc,'mmod_noinput_ps_exppar'),'png');
    
    % Pure inner dynamics without inputs and initial value -.2
    fig = tools.ParamSweep(r, [0; 3; -.2], 1, 1, 0:.1:1);
    saveFig(fig, fullfile(loc,'mmod_ps_inputpar'),'png');
    
    fig = tools.ParamSweep(r, [0; 4; -.2], 2, 1, 0:.1:1);
    saveFig(fig, fullfile(loc,'mmod_ps_inputpar2'),'png');
    
    fig = tools.ParamSweep2D(r, [0; 3; -.2], 1, [1 2], -.5:.2:1.5,-8:1.5:8);
    saveFig(fig, fullfile(loc,'mmod_ps2d_in_exp'),'png');

    function saveFig(fig, fn, fmt)
        if nargin < 3
            fmt = 'pdf';
        end
        general.Utils.saveFigure(fig, fn,'fig');
        p = get(fig,'Position');
        set(fig,'Position',[p(1:2) figsize]);
        general.Utils.saveFigure(fig, fn,fmt);
        close(fig);
    end
end