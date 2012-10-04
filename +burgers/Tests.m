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
        
        function d = getWSH12Dir
            d = fullfile(KerMor.App.DataStoreDirectory,'test_Burgers_WSH12');
        end
        
        function d = getMoRePasDir
            d = fullfile(KerMor.App.DataStoreDirectory,'test_Burgers_MoRePaS');
            if exist(d,'file') ~= 7
                [~,~] = mkdir(d);
            end
        end
        
        function createWSH12Plots(r, mu)
            % Original plot settings:
            % m.Approx.Order = [6 5];
            % r = m.buildReducedModel(44);
            
            d = models.burgers.Tests.getWSH12Dir;
            
            m = r.FullModel;
            if nargin < 2
                mu = .04;
            end
            inputidx = 1;
%             types = {'fig','pdf','jpg'};
%             types = {'png','jpg'};
            types = {'jpg'};
            
            pm = tools.PlotManager;
            pm.NoTitlesOnSave = true;
            pm.SingleSize = [720 540];
            pm.LeaveOpen = false;
            pm.UseFileTypeFolders = true;
            
            %% Simulation pics
%             m.PlotAzEl = [-49 34];
%             [~,y] = m.simulate(mu, 1);
%             [~,y1] = m.simulate(m.Data.ParamSamples(:,1), 1);
%             [~,y2] = m.simulate(m.Data.ParamSamples(:,end), 1);
%             r.ErrorEstimator.Enabled = false;
%             [t,yr] = r.simulate(mu, 1);
%             r.ErrorEstimator.Enabled = true;
%             
%             pm.FilePrefix = 'sim';
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
%             pm.savePlots(d,types,[],true);
% %             
%             %% DEIM reduction error plots
%             pm.FilePrefix = 'deim_rederr_m';
%             [errs, relerrs, times, deim_orders] = testing.DEIM.getDEIMReducedModelErrors(r, mu, inputidx);
%             testing.DEIM.getDEIMReducedModelErrors_plots(r, errs, relerrs, times, deim_orders, pm, 'estimated');
%             eold = r.ErrorEstimator;
%             r.ErrorEstimator = error.DefaultEstimator(r);
%             [terrs, trelerrs, ttimes, deim_orders] = testing.DEIM.getDEIMReducedModelErrors(r, mu, inputidx);
%             testing.DEIM.getDEIMReducedModelErrors_plots(r, terrs, trelerrs, ttimes, deim_orders, pm, 'true');
%             % Errors between both
%             testing.DEIM.getDEIMReducedModelErrors_plots(r, errs./terrs, relerrs./trelerrs, times-ttimes, deim_orders, pm, 'effectivity');
%             DEIMerrororder = r.System.f.Order(2);%#ok
%             save(fullfile(d,sprintf('reductionerrors_true_estimated_diffDEIMm_d%d',m.Dimension)),...
%                 'terrs','errs','trelerrs','relerrs','ttimes','times','deim_orders','mu','inputidx','DEIMerrororder');
%             %pm.savePlots(d,types,[1 2 5 6 9 10],true);
%             pm.savePlots(d,types);
%             pm.closeAll;
%             r.ErrorEstimator = eold;
% 
%             %% Reduction process plots
%             pm.FilePrefix = 'deim_approx';
%             h = pm.nextPlot('singvals','Singular values of SVD for DEIM basis computation');
%             semilogy(h, m.Approx.SingularValues);
%             pm.savePlots(d,types,[],true);

            %% Plot for different M,M' values
%             pm.FilePrefix = sprintf('%s_comp_m_mdash',m.SaveTag);
%             [o, eo] = meshgrid([1 2 3 4 5 7 10 12 20 25 30],...
%                 [1 2 3 4 5 7 10 15 20 25 30 40]);
%             savefile = fullfile(d,[pm.FilePrefix '.mat']);
%             if exist(savefile,'file')
%                 load(savefile);
%             else
%                 errs = testing.DEIM.computeDEIMErrors(...
%                     m.Approx, m.Data.ApproxTrainData, o(:), eo(:));
%                 save(savefile,'errs');
%             end
%             testing.DEIM.plotDEIMErrs(errs, pm);
%             figure(3); view(54,30);
%             figure(6); view(135,58);
%             pm.done;
%             pm.savePlots(d,types,[1 3 6],true);

            %% Mean required error orders over training data
            savefile = fullfile(d,sprintf('%s_minreqorders',m.SaveTag));
            if exist([savefile '.mat'],'file')
                load([savefile '.mat']);
            else
                [meanreqorder, relerrs, orders, t] = testing.DEIM.getMinRequiredErrorOrders(...
                    m, [1e-1 1e-2 1e-4 1e-6], 1:3:r.System.f.MaxOrder-1); %#ok
                t.Format = 'tex';
                t.Caption = sprintf('Minimum/mean required error orders over %d training snapshots',...
                    size(m.Data.ApproxTrainData.xi,2));
                save([savefile '.mat'],'meanreqorder','t','relerrs','orders');
            end
            t.saveToFile([savefile '.tex']);
            t.display;
%             
%             %% Estimator analysis
%             ea = tools.EstimatorAnalyzer(r);
%             ea.LineWidth = 2;
%             ea.MarkerSize = 12;
%             ea.NumMarkers = 4;
%             tmpest = r.ErrorEstimator.clone;
%             
%             %% Error estimator check for M,M' DEIM approx error est ---------------------------
%             % Using true log lip const
%             r.System.f.Order = 6;
%             tmpest.UseTrueLogLipConst = true;
%             
%             est = ea.getDefaultEstStruct;
%             est(end+1).Name = 'Reference #1'; % Expensive versions
%             est(end).Estimator = tmpest.clone;
%             est(end).Estimator.UseTrueDEIMErr = true;
%             est(end).MarkerStyle = 'p';
%             est(end).LineStyle = '-';
%             est(end).Color = [0 0.5 0];
%             ref1 = est(end);
%             est = testing.DEIM.getDEIMEstimators_ErrOrders(r,est,[1 2 3 5 10]);
%             ea.Est = est;
%             ea.SaveTexTables = fullfile(d,'table_mdash.tex');
%             [errs, relerrs, ctimes] = ea.compute(mu,1,pm);
%             ea.createPlots(errs, relerrs, ctimes, pm);
%             pm.createZoom(1,[.7 1 .9*min(errs(:,end)) 1.1*max(errs(:,end))]);
%             pm.done;
%             pm.FilePrefix = 'err_mdash_trueloglip';
%             pm.savePlots(d,types,[1 4],true);
            
%             %% Error estimator check for M,M' DEIM approx error est ---------------------------
%             % Using full jac log norm
%             r.System.f.Order = 6;
%             tmpest.UseTrueLogLipConst = false;
%             tmpest.UseFullJacobian = true;
%             ea.SaveTexTables = fullfile(d,'table_mdash.tex');
%             
%             est = ea.getDefaultEstStruct;
%             est(end+1).Name = 'Reference #2'; % Expensive versions
%             est(end).Estimator = tmpest.clone;
%             est(end).Estimator.UseTrueDEIMErr = true;
%             est(end).MarkerStyle = 'p';
%             est(end).LineStyle = '-';
%             est(end).Color = [.9 .6 0];
%             ref2 = est(end);
%             est(end+1) = ref1; % Expensive versions
%             est = testing.DEIM.getDEIMEstimators_ErrOrders(r, est, [1 2 3 5 10]);
%             ea.Est = est;
%             [errs, relerrs, ctimes] = ea.compute(mu,1,pm);
%             ea.createPlots(errs, relerrs, ctimes, pm);
%             pm.createZoom(1,[.7 1 .9*min(errs(:,end)) 1.1*max(errs(:,end))]);
%             pm.done;
%             pm.FilePrefix = 'err_mdash_fulljac_lognorm';
%             pm.savePlots(d,types,[1 4],true);
% 
%             
%             %% JacMDEIM & SimTrans error plots ------------------------------------------------
%             r.System.f.Order = [6 10];
%             tmpest.UseTrueLogLipConst = false;
%             tmpest.UseFullJacobian = false;
%             tmpest.UseTrueDEIMErr = false;
%             
%             % Expensive versions
%             est = ea.getDefaultEstStruct;
%             ref1.Estimator.UseTrueDEIMErr = false;
%             ref1.Name = 'Reference #1 m''';
%             est(end+1) = ref1;
%             ref2.Estimator.UseTrueDEIMErr = false;
%             ref2.Name = 'Reference #2 m''';
%             est(end+1) = ref2;
%             est = testing.DEIM.getDEIMEstimators_MDEIM_ST(r,est,[1 3 10],[1 3 10]);
%             ea.Est = est;
%             ea.NumMarkers = 3;
%             ea.SaveTexTables = fullfile(d,'table_jac_st.tex');
%             [~, ~, errs] = ea.start(mu,1,pm);
%             pm.createZoom(1,[.65 1 .9*min(errs(:,end)) 1.1*max(errs(:,end))]);
%             pm.done;
%             pm.FilePrefix = 'err_jac_st';
%             pm.savePlots(d,types,[1 4],true,[true false]);

%             %% Efficiencies of error estimator
%             ma = tools.ModelAnalyzer(r);
%             orders = 1:15:60;
%             for k=1:length(orders)
%                 fprintf('Using DEIM order %d...\n',orders(k));
%                 r.System.f.Order = [orders(k) 20];
%                 [errors(:,:,k), details{k}] = ma.getRedErrForParamSamples(m.Data.ParamSamples,1);%#ok
%                 [params, errors_rand(:,:,k), details_rand{k}] = ma.getRedErrForRandomParamSamples(200,1,1);%#ok
%             end
%             save(fullfile(d,sprintf('reduction_errors_d%d',m.Dimension)),...
%                 'errors','details','errors_rand','details_rand','params');
        end
        
        function createMoRePaSIIPlots(m, mu)
            % Call with d500 version
            d = models.burgers.Tests.getMoRePasDir;
            
%             types = {'fig','pdf','jpg'};
%             types = {'png','jpg'};
            types = {'jpg'};
            pm = tools.PlotManager;
            pm.NoTitlesOnSave = true;
            pm.SingleSize = [720 540];
            pm.LeaveOpen = false;
            pm.UseFileTypeFolders = true;
            
            %% Different dimension m' plots (model independent)
            load(fullfile(d,'minreqmdashorders.mat'));
            models.burgers.Tests.morepas12_plotDiffMeanReqErrs([6 7],pm,v500,v1500,v500_maxorder,v1500_maxo);
            h = legend(gca(pm.Figures(1)),'d=500','d=1500','d=500 (maxo)','d=1500 (maxo)');
            set(h,'Location','Best');
            h = legend(gca(pm.Figures(2)),'d=500','d=1500','d=500 (maxo)','d=1500 (maxo)');
            set(h,'Location','Best');
            pm.done;
            pm.FilePrefix = 'burgers';
            pm.savePlots(d,types,[],true);
            return;
            
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
%             ea = tools.EstimatorAnalyzer(r);
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
        
        function pm = wsh12_plotDiffMeanReqErrs(relerrs, orders, dims, varargin)
            % 
            % load datastore/kermor/test_Burgers_WSH12/min_mean_req_mdashorders_differentdims
            % pm = models.burgers.Tests.wsh12_plotDiffMeanReqErrs(relerrs,orders_c,dims,mo100,mo500,mo1500,mo3000);
            % 
            pm = tools.PlotManager;
            pm.LeaveOpen = true;
            nrel = size(varargin{1},2)/2;
            for k=1:length(relerrs)
                relerr = relerrs(k);
                h = pm.nextPlot(sprintf('relerr%d',k),...
                    sprintf('Minimum/Mean required m'' values for relative error < %g',relerr),...
                    'DEIM order m','m'' values');
                l = tools.LineSpecIterator(length(dims)*2);
                le = {};
                for i=1:length(dims)
                    d = varargin{i};
                    d = d(:,[k k+nrel]);
                    sel = round(d(:,1)) == d(:,1);
                    d = d(sel,:);
                    x = orders{i}(sel);
                    hold(h,'on');
                    plot(h,x,d(:,1),'Color',l.nextColor);
                    le{end+1} = sprintf('Min, d=%d',dims(i)); %#ok
                    plot(h,x,d(:,2),'Color',l.nextColor);
                    le{end+1} = sprintf('Mean, d=%d',dims(i)); %#ok
                end
                legend(h,le{:},'Location','Best');
                if k == 1
                    axis(h,[1 120 1 30]);
                end
            end
            pm.done;
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
                l = tools.LineSpecIterator(n);
                for i=1:n
                    d = varargin{i};
                    d = d(1:138,[1 idx]);
                    sel = round(d(:,1)) == d(:,1);
                    d = d(sel,:);
                    plot(h,d(:,1),d(:,2),'Color',l.nextColor,'LineWidth',2);
                end
            end
        end
        
        function m = test_Burgers
            m = models.burgers.Burgers;
            
            %% Sampling - manual
            s = sampling.ManualSampler;
            s.Samples = logspace(log10(0.04),log10(0.08),150);
            m.Sampler = s;
            
            %% Approx
            a = approx.KernelApprox;
            a.Kernel = kernels.GaussKernel;
            a.ParamKernel = kernels.GaussKernel;
            
            al = approx.algorithms.VectorialKernelOMP;
            al.MaxGFactor = [1.5 0 1.5];
            al.MinGFactor = [.2 0 .6];
            al.gameps = 1e-4;
            al.MaxExpansionSize = 300;
            al.MaxAbsErrFactor = 1e-5;
            al.MaxRelErr = 1e-3;
            al.NumGammas = 40;
            a.Algorithm = al;
            m.Approx = a;
            
            m.offlineGenerations;
            save test_Burgers;
        end
        
        function m = test_Burgers_DEIM_versions(dim, version)
            if nargin < 2
                version = 1;
                if nargin < 1
                    dim = 2000;
                end
            end
            m = models.burgers.Burgers(dim, version);
            
            if dim < 1000
                m.ModelData.TrajectoryData = data.MemoryTrajectoryData;
            else
                m.ModelData.TrajectoryData = data.FileTrajectoryData(m.ModelData);
            end
            
            %% Sampling - manual
            s = sampling.ManualSampler;
            s.Samples = logspace(log10(0.04),log10(0.08),100);
            m.Sampler = s;
            
            %% Approx
            a = approx.DEIM;
            a.MaxOrder = 80;
            m.Approx = a;
            
            offline_times = m.offlineGenerations;
            gitbranch = KerMor.getGitBranch;
            
            clear a s;
            d = fullfile(KerMor.App.DataStoreDirectory,'test_Burgers_DEIM');
            mkdir(d);
            oldd = pwd;
            cd(d);
            eval(sprintf('save test_Burgers_DEIM_d%d_v%d',dim,version));
            cd(oldd);
        end
        
        function m = getModel_WSH12(dim)
            % Original settings for experiments used in WSH12 for far.
            if nargin < 1
                dim = 100;
            end
            m = models.burgers.Burgers(dim, 2);
            %m.PlotAzEl = [-130 26];
            m.PlotAzEl = [-49 34];
            m.SaveTag = sprintf('models.burgers.d%d',m.Dimension);
            
            if dim > 500
                m.Data.TrajectoryData = [];
                m.Data.TrajectoryData = data.FileTrajectoryData(m.Data);
            end
            
            %% Sampling - log-grid
            m.System.Params(1).Range = [0.01, 0.06];
            m.System.Params(1).Desired = 100;
            s = sampling.GridSampler;
            s.Spacing = 'log';
            m.Sampler = s;
            
            %% Space reduction
            p = spacereduction.PODReducer;
            p.Mode = 'abs';
            p.Value = 44;
            m.SpaceReducer = p;
            
            %% Approx
            a = approx.DEIM;
            a.MaxOrder = 80;
            m.Approx = a;
            
            s = m.System;
            s.x0 = dscomponents.ConstInitialValue(zeros(dim,1));
            
            x = linspace(m.Omega(1),m.Omega(2),dim+2);
            x = x(2:end-1);
            pos1 = logical((x >= .1) .* (x <= .3));
            pos2 = logical((x >= .6) .* (x <= .7));
            %% Not that exciting set of inputs
%             s.Inputs{1} = @(t)2*sin(2*t*pi);
%             B = zeros(dim,1);
%             B(pos1) = 4*exp(-((x(pos1)-.2)/.03).^2);
%             B(pos2) = 1;
%             
            s.Inputs{1} = @(t)[sin(2*t*pi); (t>.2)*(t<.4)];
            B = zeros(dim,2);
            B(pos1,1) = 4*exp(-((x(pos1)-.2)/.03).^2);
            B(pos2,2) = 4;
            
            s.B = dscomponents.LinearInputConv(B);

            m.TrainingInputs = 1;
            
%             t(1) = m.off1_createParamSamples;
%             t(2) = m.off2_genTrainingData;
%             t(3) = m.off3_computeReducedSpace;
%             t(4) = m.off4_genApproximationTrainData;
%             t(5) = m.off5_computeApproximation;
%             t(6) = m.off6_prepareErrorEstimator;
%             offline_times = t;
            offline_times = m.offlineGenerations;
            gitcommit = KerMor.getGitBranch;
            
            clear a s;
            d = models.burgers.Tests.getWSH12Dir;
            [~,~] = mkdir(d);
            oldd = pwd;
            cd(d);
            eval(sprintf('save test_Burgers_WSH12_d%d',dim));
            cd(oldd);
        end
        
        function m = getModel_WSH12_withextras(dim, withfxi, withfd, withbspan)
            if nargin < 4
                withbspan = false;
                if nargin < 3
                    withfd = false;
                    if nargin < 2
                        withfxi = false;
                        if nargin < 1
                            dim = 100;
                        end
                    end
                end
            end
            m = models.burgers.Burgers(dim, 2);
            %m.PlotAzEl = [-130 26];
            m.PlotAzEl = [-49 34];
            m.SaveTag = sprintf('models.burgers.d%d.fx%d.fd%d.bs%d',...
                m.Dimension,withfxi,withfd,withbspan);
            
            if dim > 500
                m.Data.TrajectoryData = [];
                m.Data.TrajectoryData = data.FileTrajectoryData(m.Data);
            end
            
            %% Sampling - log-grid
            m.System.Params(1).Range = [0.01, 0.06];
            m.System.Params(1).Desired = 100;
            s = sampling.GridSampler;
            s.Spacing = 'log';
            m.Sampler = s;
            
            %% Space reduction
            m.ComputeTrajectoryFxiData = withfxi;
            m.ComputeBSpan = withbspan;
            p = spacereduction.PODReducer;
            p.IncludeTrajectoryFxiData = withfxi;
            p.IncludeFiniteDifferences = withfd;
            p.IncludeBSpan = withbspan;
            p.Mode = 'abs';
            p.Value = 200;
            m.SpaceReducer = p;
            
            %% Approx
            a = approx.DEIM;
            a.MaxOrder = 200;
            m.Approx = a;
            
            s = m.System;
            s.x0 = dscomponents.ConstInitialValue(zeros(dim,1));
            
            x = linspace(m.Omega(1),m.Omega(2),dim+2);
            x = x(2:end-1);
            pos1 = logical((x >= .1) .* (x <= .3));
            pos2 = logical((x >= .6) .* (x <= .7));
            %% Not that exciting set of inputs
%             s.Inputs{1} = @(t)2*sin(2*t*pi);
%             B = zeros(dim,1);
%             B(pos1) = 4*exp(-((x(pos1)-.2)/.03).^2);
%             B(pos2) = 1;
%             
            s.Inputs{1} = @(t)[sin(2*t*pi); (t>.2)*(t<.4)];
            B = zeros(dim,2);
            B(pos1,1) = 4*exp(-((x(pos1)-.2)/.03).^2);
            B(pos2,2) = 4;
            
            s.B = dscomponents.LinearInputConv(B);

            m.TrainingInputs = 1;
%             t(1) = m.off1_createParamSamples;
%             t(2) = m.off2_genTrainingData;
%             t(3) = m.off3_computeReducedSpace;
%             t(4) = m.off4_genApproximationTrainData;
%             t(5) = m.off5_computeApproximation;
%             t(6) = m.off6_prepareErrorEstimator;
%             offline_times = t;
            offline_times = m.offlineGenerations;
            gitcommit = KerMor.getGitBranch;
            
            clear a s;
            d = models.burgers.Tests.getWSH12Dir;
            [~,~] = mkdir(d);
            oldd = pwd;
            cd(d);
            extra = '';
            if withfxi
                extra = '_withFXI';
            end
            if withfd
                extra = [extra '_withFD'];
            end
            if withbspan
                extra = [extra '_withBS'];
            end
            eval(sprintf('save test_Burgers_WSH12_d%d%s',dim,extra));
            cd(oldd);
        end
    end
    
end