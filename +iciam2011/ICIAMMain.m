function [m, r, d] = ICIAMMain(dim)
    % Current version works with KerMor 0.5

%     clear classes;
    if nargin == 0
        dim = 1000;
    end
    %%
    mu = [.5; 1];
    in = 1;
    
    m = models.iciam2011.ICIAMExperiment(dim);
    m.offlineGenerations;
    r = m.buildReducedModel;
    
    %% Analyze
    d = EstimatorAnalyzer;
    d.EstimatorIterations = [1 2];
    d.EstimatorVersions = [1 1 0 0 1 0 0 1 1];
    %d.SingleFigures = true;
    d.setModel(r);
%     
    d.compute(mu,in);
%     
%     ma = ModelAnalyzer;
%     ma.SingleFigures = d.SingleFigures;
%     ma.analyze(r, mu, in);
    
%     ParamSweep(r, mu, in, 1, 0:.1:.2); %-.9:.5:1.5
    [pm,Y,E] = ParamSweep(r, mu, in, 1, 0:.1:1);
    
    %% Experiments regardign the initial value
    [pm,Y,E] = ParamSweep(r, [1; 1], in, 2, 0.25:.1:1.2);
    [pm,Y,E] = ParamSweep(r, [.5; 1], in, 2, 0.25:.1:1.2);
    [pm,Y,E] = ParamSweep(r, [0; 1], in, 2, 0.25:.1:1.2);
    % Extranah
%     ParamSweep(r, [.7; 1], 1, 2, 0.9:.001:.94);
end