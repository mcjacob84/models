classdef Tests
% Tests: Some test settings regarding the PCD model simulations.
%
%
%
% @author Daniel Wirtz @date 2012-06-11
%
% @new{0,6,dw,2012-06-11} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods(Static)
        
        %% ---------------- 2D tests --------------------
        function reductionErrorAnalysis_2D(d, r, pm)
            % Reduction error analysis - computes the reduction errors for 
            % training parameters and random parameters

            m = r.FullModel;
            ma = tools.ModelAnalyzer(r);
            %orders = [3 15 35 55 100]; % orders for old large 3000s model
            orders = [1   5   20  97 102 107 145;...
                      149 145 130 53 48  43  5];
            errors = zeros(8,m.Data.SampleCount,length(orders));
            ns = 100;
            rand_errors = zeros(8,ns,length(orders));
            details = cell(1,length(orders));
            rand_details = details;
            for j=1:size(orders,2)
                fprintf('Setting DEIM order = %d\n',orders(1,j));
                r.System.f.Order = orders(:,j)';
                [e, details{j}] = ma.getRedErrForParamSamples;
                errors(:,:,j) = e;
                [params, re, rand_details{j}] = ma.getRedErrForRandomParamSamples(ns,1); % 100 params, seed 1
                rand_errors(:,:,j) = re;
            end
            save(fullfile(d,'reduction_errors.mat'),'errors','orders','details');
            save(fullfile(d,'reduction_errors_randparams.mat'),'rand_errors','orders','params','rand_details');
            
%             s=load(fullfile(d,'reduction_errors.mat'));
%             sr=load(fullfile(d,'reduction_errors_randparams.mat'));
%             sr.rand_errors = sr.errors; % fix for re-runs
%             orders = s.orders;
%             
%             str = sprintf(' on %d training parameters, m=%d',size(s.errors,2),sr.orders(end));
%             strr = sprintf(' on %d random parameters, m=%d',size(sr.rand_errors,2),sr.orders(end));
%             directplot(s.errors,'train',str);
%             directplot(sr.rand_errors,'random',strr);
%             paramplot(s.errors,'train',str,m.Data.ParamSamples,7);
%             paramplot(sr.rand_errors,'random',strr,sr.params,7);
%             pm.done;
            
%             function directplot(errors, tag, str)
%                 h = pm.nextPlot([tag '_abs'],['Absolute Linf-L2 reduction errors' str],...
%                 'DEIM order m','param sample nr');
%                 tools.LogPlot.logsurf(h, orders, 1:size(errors,2), squeeze(errors(1,:,:)));
%                 h = pm.nextPlot([tag '_rel'],['Relative Linf-L2 reduction errors' str],...
%                     'DEIM order m','param sample nr');
%                 tools.LogPlot.logsurf(h, orders, 1:size(errors,2), squeeze(errors(2,:,:)));
%             end
%             
%             function paramplot(errors, tag, str, params, orderidx)
%                 innerPlot(1,'abs');
%                 innerPlot(2,'rel');
%                 function innerPlot(pos,tag2)
%                     h = pm.nextPlot(['paramdomain_' tag '_' tag2],...
%                         ['Locations of 10% of worst ' tag2 ' errors' str],...
%                     ['param "' m.System.Params(1).Name '"'],...
%                     ['param "' m.System.Params(2).Name '"']);
%                     
%                     [v, idx] = sort(errors(pos,:,orderidx));
%                     v = log10(v);
%                     merr = min(v);
%                     range = (max(v)-merr);
%                     hold(h,'on');
%                     for k=1:length(v)
%                         if k > length(v)*.85
%                             c2 = 'blue';
%                         else
%                             c2 = 'none';
%                         end
%                         c = [(v(k)-merr)/range 1-(v(k)-merr)/range 0];
%                         plot3(h,params(1,idx(k)),params(2,idx(k)),v(k),...
%                             'MarkerFaceColor',c,...%'Color',c,...
%                             'MarkerEdgeColor',c2,...
%                             'Marker','o',...
%                             'MarkerSize',10);
%                     end
%                     view(h,21,24);
%                     hold(h,'off');
%                 end
%             end
        end
      
        %% Model creation functions
        function m = tests_PCD_DEIM_2D(dim)
            % Original setting for the large-scale reduction up to T=3000 with 200 lin-spaced
            % parameters and 120'DEIM, 80'JacMDEIM
            %
            % Computed partial similarity transform with target size 50 over 20000
            % eigenvectors. Resulting reduction 50/60000 (99.9167%) 
            m = models.pcd.PCDModel(2);
            
            m.T = 3000; %[s]
            m.dt = 5; %[s]
            if nargin < 1
                dim = 150;
            end
            m.System.h = (m.System.Omega(1,2)-m.System.Omega(1,1))/(dim-1);
            
            if all([m.System.f.fDim m.System.f.xDim] < 1000)
                m.Data.TrajectoryData = data.MemoryTrajectoryData;
            end
            
            % area
            m.System.Params(1).Desired = 10;
            % rate
            m.System.Params(2).Desired = 20;
            s = sampling.GridSampler;
            s.Spacing = 'lin';
            m.Sampler = s;
            
            s = spacereduction.PODGreedy;
            s.Eps = 1e-7;
            m.SpaceReducer = s;
            
            m.Approx = approx.DEIM;
            m.Approx.MaxOrder = 120;
            
            s = data.selection.LinspaceSelector;
            s.Size = 20000;
            m.Approx.TrainDataSelector = s;
            
            m.System.MaxTimestep = m.dt;
            m.ODESolver = solvers.SemiImplicitEuler(m);
            
            e = error.DEIMEstimator;
            e.JacMatDEIMMaxOrder = 80;
            ts = data.selection.LinspaceSelector;
            ts.Size = round(.2 * m.T * m.System.Params(1).Desired / m.dt);
            e.TrainDataSelector = ts;
            m.ErrorEstimator = e;
            
            gitbranch = KerMor.getGitBranch;
            d = fullfile(KerMor.App.DataStoreDirectory,sprintf('tests_PCD_DEIM_2D_%d_%d',m.System.Dims));
            [~,~] = mkdir(d);
            oldd = pwd;
            cd(d);
            save tests_PCD_DEIM_2D;
            cd(oldd);
            
%             t(1) = m.off1_createParamSamples;
%             t(2) = m.off2_genTrainingData;
%             t(3) = m.off3_computeReducedSpace;
%             t(4) = m.off4_genApproximationTrainData;
%             t(5) = m.off5_computeApproximation;
%             t(6) = m.off6_prepareErrorEstimator;
%             offline_times = t;

            offline_times = m.offlineGenerations;
            
            clear s t;
            oldd = pwd;
            cd(d);
            save tests_PCD_DEIM_2D;
            cd(oldd);
        end
        
        function m = tests_PCD_DEIM_2D_500s(dim)
            % New configuration with shorter runtime
            % Samples are randomly distributed!
            m = models.pcd.PCDModel(2);
            
            m.T = 500; %[s]
            m.dt = 5; %[s]
            if nargin < 1
                dim = 150;
            end
            m.System.h = (m.System.Omega(1,2)-m.System.Omega(1,1))/(dim-1);
            
%             if all([m.System.f.fDim m.System.f.xDim] < 1000)
%                 m.Data.TrajectoryData = data.MemoryTrajectoryData;
%             end
            
            s = sampling.RandomSampler;
            s.Seed = 35178;
            s.Samples = 200;
            m.Sampler = s;
            
            m.ComputeTrajectoryFxiData = true;
            s = spacereduction.PODGreedy;
            s.Eps = 1e-7;
            s.MaxSubspaceSize = 200;
            s.MinRelImprovement = 1e-8;
            s.IncludeTrajectoryFxiData = true; % enable inclusion of f evaluations!
            m.SpaceReducer = s;
            
            m.Approx = approx.DEIM;
            m.Approx.MaxOrder = 250;
            
            s = data.selection.LinspaceSelector;
            s.Size = 30000;
            m.Approx.TrainDataSelector = s;
            
            m.System.MaxTimestep = m.dt;
            m.ODESolver = solvers.SemiImplicitEuler(m);
            
            e = error.DEIMEstimator;
            e.JacMatDEIMMaxOrder = 200;
            e.JacSimTransMaxSize = 200;
            ts = data.selection.LinspaceSelector;
            ts.Size = 30000;
            e.TrainDataSelector = ts;
            m.ErrorEstimator = e;
            
            gitbranch = KerMor.getGitBranch;
            d = fullfile(KerMor.App.DataStoreDirectory,sprintf('tests_PCD_DEIM_2D_%d_%d_500s',m.System.Dims));
            [~,~] = mkdir(d);
            oldd = pwd;
            cd(d);
            save tests_PCD_DEIM_2D_500s;
            cd(oldd);
            
%             t(1) = m.off1_createParamSamples;
%             t(2) = m.off2_genTrainingData;
%             t(3) = m.off3_computeReducedSpace;
%             t(4) = m.off4_genApproximationTrainData;
%             t(5) = m.off5_computeApproximation;
%             t(6) = m.off6_prepareErrorEstimator;
%             offline_times = t;

            offline_times = m.offlineGenerations;
            
            clear s t;
            oldd = pwd;
            cd(d);
            save tests_PCD_DEIM_2D_500s;
            cd(oldd);
        end
        
        function m = tests_PCD_DEIM_1D(dim)
            m = models.pcd.PCDModel(1);
            
            m.T = 3000; %[s]
            m.dt = 5; %[s]
            if nargin < 1
                dim = 25;
            end
            m.System.h = (m.System.Omega(2)-m.System.Omega(1))/(dim-1);
            
            m.Data.TrajectoryData = data.MemoryTrajectoryData;
            
            m.System.Params(1).Desired = 20;
            s = sampling.GridSampler;
            s.Spacing = 'lin';
            m.Sampler = s;
            
            m.Approx = approx.DEIM;
            m.Approx.MaxOrder = 60;
            
            s = data.selection.DefaultSelector;
            %s.Size = 15000;
            m.Approx.TrainDataSelector = s;
            
            m.System.MaxTimestep = m.dt;
            m.ODESolver = solvers.SemiImplicitEuler(m);
            
            e = error.DEIMEstimator;
            e.JacMatDEIMMaxOrder = 60;
            ts = data.selection.LinspaceSelector;
            ts.Size = round(.2 * m.T * m.System.Params(1).Desired / m.dt);
            e.TrainDataSelector = ts;
            m.ErrorEstimator = e;
            
            offline_times = m.offlineGenerations;
            gitbranch = KerMor.getGitBranch;
            
            clear s;
            d = fullfile(KerMor.App.DataStoreDirectory,'tests_PCD_DEIM_1D');
            [~,~] = mkdir(d);
            oldd = pwd;
            cd(d);
            %eval(sprintf('save tests_PCD_DEIM_1D',dim,version));
            save tests_PCD_DEIM_1D;
            cd(oldd);
        end 
    end   
end