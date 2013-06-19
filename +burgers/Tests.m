classdef Tests
% Tests: Some tests and simulation settings for Burger's equation models.
%
%
%
% @author Daniel Wirtz @date 2012-06-13
%
% @new{0,6,dw,2012-06-13} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods(Static)
        
        function d = getMoRePasDir
            d = fullfile(KerMor.App.DataDirectory,'test_Burgers_MoRePaS');
            if exist(d,'file') ~= 7
                [~,~] = mkdir(d);
            end
        end
        
        function createMoRePaSIIPlots(m, mu)
            % Call with d500 version
            %
            % Parameters:
            % m: The model @type models.Burgers.Burgers
            d = models.burgers.Tests.getMoRePasDir;
            
%             types = {'fig','pdf','jpg'};
%             types = {'png','jpg'};
            types = {'jpg'};
            pm = PlotManager;
            pm.NoTitlesOnSave = true;
            pm.SingleSize = [720 540];
            pm.LeaveOpen = false;
            pm.UseFileTypeFolders = true;
            
            %% Different dimension m' plots (model independent)
%             load(fullfile(d,'minreqmdashorders.mat'));
%             models.burgers.Tests.morepas12_plotDiffMeanReqErrs([6 7],pm,v500,v1500,v500_maxorder,v1500_maxo);
%             h = legend(gca(pm.Figures(1)),'d=500','d=1500','d=500 (maxo)','d=1500 (maxo)');
%             set(h,'Location','Best');
%             h = legend(gca(pm.Figures(2)),'d=500','d=1500','d=500 (maxo)','d=1500 (maxo)');
%             set(h,'Location','Best');
%             pm.done;
%             pm.FilePrefix = 'burgers';
%             pm.savePlots(d,types,[],true);
%             return;
            
            %% Model preparations
            %m = r.FullModel;
            m.Approx.Order = [12 5];
            r = m.buildReducedModel(100);
            
            if nargin < 2
                mu = .0400574;
            end
            inputidx = 1;
            
            %% Simulation pics
%             m.PlotAzEl = [-49 34];
%             [~,y] = m.simulate(mu, 1);
%             [~,y1] = m.simulate(m.Data.ParamSamples(:,1), 1);
%             [~,y2] = m.simulate(m.Data.ParamSamples(:,end), 1);
%             r.ErrorEstimator.Enabled = false;
%             [t,yr] = r.simulate(mu, 1);
%             r.ErrorEstimator.Enabled = true;
%             
%             pm.FilePrefix = [m.SaveTag '_sim'];
%             m.plot(t, y1, pm.nextPlot('full_sol_mumin',...
%                 sprintf('Full solution for \\mu=%g',m.Data.ParamSamples(:,1))));
%             
%             m.plot(t, y2, pm.nextPlot('full_sol_mumax',...
%                 sprintf('Full solution for \\mu=%g',m.Data.ParamSamples(:,end))));
%             %view(-24,44);
%             m.plot(t, y1-y2, pm.nextPlot('full_sol_muvar',...
%                 sprintf('Model variance over parameter range [%g, %g]',...
%                 m.Data.ParamSamples(:,1),m.Data.ParamSamples(:,end))));
%             %view(-26,22);
%             m.plot(t, y, pm.nextPlot('full_sol',...
%                 sprintf('Full solution for \\mu=%g',mu)));
%             m.plot(t, y, pm.nextPlot('red_sol',...
%                 sprintf('Reduced solution for \\mu=%g',mu)));
%             m.plot(t, y-yr, pm.nextPlot('abs_err','Absolute error'));
% %             m.plot(t, abs((y-yr)./y), pm.nextPlot('rel_err','Relative error'));
%             pm.done;
%             pm.savePlots(d,types,[]);

            % Plot for different M,M' values
%             pm.FilePrefix = sprintf('%s_m_mdash_errors',m.SaveTag);
%             savefile = fullfile(d,pm.FilePrefix);
%             if exist([savefile '.mat'],'file')
%                 load([savefile '.mat']);
%             else
%                 errordata = testing.DEIM.computeDEIMErrors(...
%                     m.Approx, m.Data.ApproxTrainData);
%                 relerrs = [1e-1 1e-2 1e-4 1e-6];
%                 [t, valtrue] = testing.DEIM.getMinReqErrorOrdersTable(errordata, relerrs, m.Data.ApproxTrainData.xi.m);
%                 t.saveToFile([savefile '.tex']);
%                 [tmaxo, valmaxorder] = testing.DEIM.getMinReqErrorOrdersTable(errordata, relerrs, m.Data.ApproxTrainData.xi.m, m.Approx.MaxOrder);
%                 tmaxo.saveToFile([savefile '_true.tex']);
%                 save([savefile '.mat'],'errordata', 'relerrs', 't', 'tmaxo','valtrue','valmaxorder');
%             end
%             
%             t.display;
%             tmaxo.display;
%             testing.DEIM.plotDEIMErrs(errordata, pm);
%             figure(2); view(28,12);
%             figure(6); view(-34,46);
%             pm.done;
%             pm.savePlots(d,types,[1 2 6],true);

            %% Estimator analysis
%             ea = EstimatorAnalyzer(r);
%             ea.LineWidth = 2;
%             ea.MarkerSize = 12;
%             ea.NumMarkers = 4;
%             ea.PlotStartIndex = 4;
%             e = r.ErrorEstimator;
            
            %% Error estimator check for M,M' DEIM approx error est ---------------------------
            % Using true log lip const
%             e.UseTrueLogLipConst = true;
%             
%             est = ea.getDefaultEstStruct;
%             est(end+1).Name = 'Ref: True local const'; % Expensive versions
%             est(end).Estimator = e.clone;
%             est(end).Estimator.UseTrueDEIMErr = true;
%             est(end).MarkerStyle = 'p';
%             est(end).LineStyle = '-';
%             est(end).Color = [0 0.5 0];
%             ref1 = est(end);
%             est = testing.DEIM.getDEIMEstimators_ErrOrders(r,est,[1 2 3 5 10]);
%             ea.Est = est;
%             ea.SaveTexTables = fullfile(d,'table_mdash.tex');
%             [errs, relerrs, ctimes] = ea.compute(mu,1);
%             ea.createPlots(errs, relerrs, ctimes, pm);
%             pm.createZoom(1,[.7 1 .9*min(errs(:,end)) 1.1*max(errs(:,end))]);
%             pm.done;
%             pm.FilePrefix = sprintf('%s_err_mdash_trueloglip',m.SaveTag);
%             pm.savePlots(d,types,[1 4],true);
            
%             %% Error estimator check for M,M' DEIM approx error est ---------------------------
%             % Using full jac log norm
%             e.UseTrueLogLipConst = false;
%             e.UseFullJacobian = true;
%             ea.SaveTexTables = fullfile(d,'table_mdash.tex');
%             
%             est = ea.getDefaultEstStruct;
%             est(end+1).Name = 'L_G[J(y)]'; % Expensive versions
%             est(end).Estimator = e.clone;
%             est(end).Estimator.UseTrueDEIMErr = true;
%             est(end).MarkerStyle = 'p';
%             est(end).LineStyle = '-';
%             est(end).Color = [.9 .6 0];
%             ref2 = est(end);
%             est(end+1) = ref1; % Expensive versions
%             %est = testing.DEIM.getDEIMEstimators_ErrOrders(r, est, [1 2 3 5 10]);
%             ea.Est = est;
%             [errs, relerrs, ctimes] = ea.compute(mu,1);
%             ea.createPlots(errs, relerrs, ctimes, pm);
%             pm.createZoom(1,[.7 1 .9*min(errs(:,end)) 1.1*max(errs(:,end))]);
%             pm.done;
%             pm.FilePrefix = sprintf('%s_err_mdash_fulljac_lognorm',m.SaveTag);
%             pm.savePlots(d,types,[1 4],true);
% 
%             
%             %% JacMDEIM & SimTrans error plots ------------------------------------------------
%             r.System.f.Order = [12 10];
%             e.UseTrueLogLipConst = false;
%             e.UseFullJacobian = false;
%             e.UseTrueDEIMErr = false;
%             
%             % Expensive versions
%             est = ea.getDefaultEstStruct;
%             ref1.Estimator.UseTrueDEIMErr = false;
%             %ref1.Name = 'Reference #1 m''';
%             est(end+1) = ref1;
%             ref2.Estimator.UseTrueDEIMErr = false;
%             %ref2.Name = 'Reference #2 m''';
%             est(end+1) = ref2;
%             est = testing.DEIM.getDEIMEstimators_MDEIM_ST(r,est,[1 3 10],[1 3 10]);
%             ea.Est = est;
%             ea.NumMarkers = 3;
%             ea.SaveTexTables = fullfile(d,'table_jac_st.tex');
%             [errs, relerrs, ctimes] = ea.compute(mu,1);
%             ea.createPlots(errs, relerrs, ctimes, pm);
%             pm.createZoom(1,[.65 1 .9*min(errs(:,end)) 1.1*max(errs(:,end))]);
%             pm.done;
%             pm.FilePrefix = 'err_jac_st';
%             pm.savePlots(d,types,[1 4],true,[true false]);
        end
        
        function morepas12_plotDiffMeanReqErrs(indices, pm, varargin)
            % 
            % load minreqmdashorders
            % pm = models.burgers.Tests.morepas12_plotDiffMeanReqErrs(6,[100 500 101 501],v100,v500,v100_maxorder,v500_maxorder);
            % 
            for k = 1:length(indices)
                idx = indices(k);
                h = pm.nextPlot(sprintf('diffmeanreqerrs_%d',idx),...
                        'Minimum/Mean required m'' values for relative errors',...
                        'DEIM order m','m'' values');
                hold(h,'on');
                n = length(varargin);
                l = LineSpecIterator(n);
                for i=1:n
                    d = varargin{i};
                    d = d(1:138,[1 idx]);
                    sel = round(d(:,1)) == d(:,1);
                    d = d(sel,:);
                    plot(h,d(:,1),d(:,2),'Color',l.nextColor,'LineWidth',2);
                end
            end
        end
        
        function [res, m] = test_Burgers_DEIM_versions(dim, version)
            if nargin < 2
                version = 1;
                if nargin < 1
                    dim = 20;
                end
            end
            m = models.burgers.Burgers(dim, version);
            
            if dim > 1000
                m.Data.useFileTrajectoryData;
            end
            
            %% Sampling - manual
            s = sampling.ManualSampler;
            s.Samples = logspace(log10(0.04),log10(0.08),10);
            m.Sampler = s;
            
            %% Approx
            a = approx.DEIM;
            a.MaxOrder = 40;
            m.Approx = a;
            
            m.offlineGenerations;
            a.Order = [5 2];
            r = m.buildReducedModel;
            [t,y] = r.simulate(r.getRandomParam);
            m.plot(t,y);
            res = true;
        end
    end
    
end