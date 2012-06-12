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
        %% ---------------- 1D tests --------------------
        
        function m = tests_PCD_DEIM_1D
            m = models.pcd.PCDModel(1);
            
            m.T = 10000; %[s]
            m.dt = 5; %[s]
            
            m.Data = data.MemoryModelData;
            
            m.System.Params(1).Desired = 40;
            s = sampling.GridSampler;
            s.Spacing = 'lin';
            m.Sampler = s;
            
            m.Approx = approx.DEIM;
            m.Approx.MaxOrder = 60;
            
            s = approx.selection.DefaultSelector;
            %s.Size = 15000;
            m.Approx.TrainDataSelector = s;
            
            m.System.MaxTimestep = m.dt;
            m.ODESolver = solvers.ode.SemiImplicitEuler(m);
            
            e = error.DEIMEstimator;
            e.JacMatDEIMMaxOrder = 60;
            ts = approx.selection.LinspaceSelector;
            ts.Size = round(.2 * m.T * m.System.Params(1).Desired / m.dt);
            e.TrainDataSelector = ts;
            m.ErrorEstimator = e;
            
            offline_times = m.offlineGenerations;
            gitbranch = KerMor.getGitBranch;
            
            clear s;
            d = fullfile(KerMor.App.DataStoreDirectory,'tests_PCD_DEIM_1D');
            mkdir(d);
            oldd = pwd;
            cd(d);
            %eval(sprintf('save tests_PCD_DEIM_1D',dim,version));
            save tests_PCD_DEIM_1D;
            cd(oldd);
        end
        
        function m = testConfig1
            % 1D pcd model, approx possibility test 1
            %
            % - T = 12
            % - dt = .01
            % - memory model data
            % - linspace selector with 20000 pts
            % - linear parameter space sampling.
            % - No subspace projection.
            
            % - default selector (=all traj data for learning),
            m = models.pcd.PCDModel(1);
            
            m.T = 12; %[s]
            m.dt = .01; %[s]
            
            s = sampling.GridSampler;
            s.Spacing = 'lin';
            m.Sampler = s;
            
            a = m.Approx;
            %s = approx.selection.DefaultSelector;
            s = approx.selection.LinspaceSelector;
            s.Size = 20000;
            a.TrainDataSelector = s;
            aa = approx.algorithms.AdaptiveCompWiseKernelApprox;
            aa.MaxExpansionSize = 300;
            aa.MaxRelErr = 1e-5;
            aa.MaxAbsErrFactor = 1e-5;
            aa.MaxGFactor = 5;
            aa.MinGFactor = 0.001;
            aa.NumGammas = 60;
            a.Algorithm = aa;
            
            m.Data = data.MemoryModelData;
            
            m.SpaceReducer = [];
        end
        
        function m = testConfig2
            % Same as testConfig1, but with zero initial conditions.
            
            m = models.pcd.PCDModel.testConfig1;
            
            dim = size(m.System.x0.evaluate([]),1);
            m.System.x0 = dscomponents.ConstInitialValue(zeros(dim,1));
        end
        
        function oldGoodConfigSearch
            % Configuration to try and find the good old approximation with
            % little approximation error
            m = models.pcd.PCDModel(1);
            
            m.T = 6; %[s]
            m.dt = .001; %[s]
            m.System.h = m.L/24;
            o = solvers.ode.MLode15i;
            o.AbsTol = 1e-6;
            o.RelTol = 1e-5;
            o.MaxStep = [];
            m.ODESolver = o;
            
            m.System.Params(1).Range = [1e-3 .1];
            m.System.Params(1).Desired = 12;
            
            s = sampling.GridSampler;
            s.Spacing = 'lin';
            m.Sampler = s;
            
            a = m.Approx;
%             s = approx.selection.DefaultSelector;
%             s = approx.selection.LinspaceSelector;
            s = approx.selection.TimeSelector;
            s.Size = 12000;
            a.TrainDataSelector = s;
            aa = approx.algorithms.AdaptiveCompWiseKernelApprox;
            aa.MaxExpansionSize = 200;
            aa.MaxRelErr = 1e-5;
            aa.MaxAbsErrFactor = 1e-5;
            a.Algorithm = aa;
            
            m.Data = data.MemoryModelData;
            
            % Zero initial conditions
            dim = size(m.System.x0.evaluate([]),1);
            m.System.x0 = dscomponents.ConstInitialValue(zeros(dim,1));
            
            m.SpaceReducer = [];
            
            t1 = m.off1_createParamSamples;
            t2 = m.off2_genTrainingData;
            t3 = m.off3_computeReducedSpace;
            t4 = m.off4_genApproximationTrainData;
            
            factors = general.Utils.createCombinations([10 15 20], [1 2 3 4 5],[.1 .01 .001]);
            n = size(factors,2);
            K = approx.KernelApprox.empty;
            t5 = zeros(1,n);
            for i = 1:n
                fprintf('Starting approximation run %d of %d..\n',i,n);
                aa.NumGammas = factors(1,i);
                aa.MaxGFactor = factors(2,i);
                aa.MinGFactor = factors(3,i);
                t5(i) = m.off5_computeApproximation;
                K(i) = m.Approx.clone;
            end
            times = [t1 t2 t3 t4 t5];
            
            save oldConfigSearch m K times factors;
        end
        
        function oldGoodConfigSearch_Large
            % Configuration to try and find the good old approximation with
            % little approximation error
            m = models.pcd.PCDModel(1);
            
            m.T = 6; %[s]
            m.dt = .001; %[s]
            m.System.h = m.L/24;
            o = solvers.ode.MLode15i;
            o.AbsTol = 1e-6;
            o.RelTol = 1e-5;
            o.MaxStep = [];
            m.ODESolver = o;
            
            m.System.Params(1).Range = [1e-3 .1];
            m.System.Params(1).Desired = 20;
            
            s = sampling.GridSampler;
            s.Spacing = 'lin';
            m.Sampler = s;
            
            a = m.Approx;
%             s = approx.selection.DefaultSelector;
%             s = approx.selection.LinspaceSelector;
            s = approx.selection.TimeSelector;
            s.Size = 24000;
            a.TrainDataSelector = s;
            aa = approx.algorithms.AdaptiveCompWiseKernelApprox;
            aa.MaxExpansionSize = 300;
            aa.MaxRelErr = 1e-5;
            aa.MaxAbsErrFactor = 1e-5;
            a.Algorithm = aa;
            
            m.Data = data.MemoryModelData;
            
            % Zero initial conditions
            dim = size(m.System.x0.evaluate([]),1);
            m.System.x0 = dscomponents.ConstInitialValue(zeros(dim,1));
            
            m.SpaceReducer = [];
            
            t1 = m.off1_createParamSamples;
            t2 = m.off2_genTrainingData;
            t3 = m.off3_computeReducedSpace;
            t4 = m.off4_genApproximationTrainData;
            
            factors = general.Utils.createCombinations([5 10 15 20], [1 2 3 4 5],[.1 .01 .001]);
            n = size(factors,2);
            K = approx.KernelApprox.empty;
            t5 = zeros(1,n);
            for i = 1:n
                fprintf('Starting approximation run %d of %d..\n',i,n);
                aa.NumGammas = factors(1,i);
                aa.MaxGFactor = factors(2,i);
                aa.MinGFactor = factors(3,i);
                t5(i) = m.off5_computeApproximation;
                K(i) = m.Approx.clone;
            end
            times = [t1 t2 t3 t4 t5];
            
            save oldConfigSearch_Large m K times factors;
        end
        
        %% ----------- Kernel approximation tests -----------
        function r = oldconf_search_VKOGA
            m = models.pcd.PCDModel(1);
            m.T = 6;
            m.dt = .001;
            m.System.Params(1).Range = [1e-3 .1];
            m.System.Params(1).Desired = 12;
            
            s = sampling.GridSampler;
            s.Spacing = 'lin';
            m.Sampler = s;
            
            a = approx.KernelApprox;
            a.TimeKernel = kernels.NoKernel;%kernels.GaussKernel(1);
            %a.TimeKernel.G = 1;
            a.Kernel = kernels.GaussKernel(1);
            a.Kernel.G = 1;
            a.ParamKernel = kernels.GaussKernel(1);
            a.ParamKernel.G = 1;
            
            s = approx.selection.LinspaceSelector;
            s.Size = 12000;
            a.TrainDataSelector = s;
            
            aa = approx.algorithms.VectorialKernelOMP;
            aa.MaxExpansionSize = 200;
            aa.MaxRelErr = 1e-5;
            aa.MaxAbsErrFactor = 1e-5;
            aa.MaxGFactor = 50;
            aa.MinGFactor = .01;
            aa.NumGammas = 40;
            a.Algorithm = aa; 
            
            m.Approx = a;
            
            m.offlineGenerations;
            r = m.buildReducedModel;
        end
        
        function nonzerox0_300s_dt01_largeAA_subsp
            % Configuration to try and find the good old approximation with
            % little approximation error
            m = models.pcd.PCDModel(1);
            
            m.T = 300; %[s]
            m.dt = .1; %[s]
            m.System.h = m.L/24;
            o = solvers.ode.MLode15i;
            o.AbsTol = 1e-6;
            o.RelTol = 1e-5;
            o.MaxStep = [];
            m.ODESolver = o;
            
            m.System.Params(1).Range = [1e-3 .1];
            m.System.Params(1).Desired = 20;
            
            s = sampling.GridSampler;
            s.Spacing = 'lin';
            m.Sampler = s;
            
            a = m.Approx;
%             s = approx.selection.DefaultSelector;
%             s = approx.selection.LinspaceSelector;
            s = approx.selection.TimeSelector;
            s.Size = 24000;
            a.TrainDataSelector = s;
            aa = approx.algorithms.AdaptiveCompWiseKernelApprox;
            aa.MaxExpansionSize = 300;
            aa.MaxRelErr = 1e-5;
            aa.MaxAbsErrFactor = 1e-5;
            a.Algorithm = aa;
            
            m.Data = data.MemoryModelData;
            
            % Zero initial conditions
            %dim = size(m.System.x0.evaluate([]),1);
            %m.System.x0 = dscomponents.ConstInitialValue(zeros(dim,1));
            
            %m.SpaceReducer = [];
            m.SpaceReducer = spacereduction.PODGreedy;
            m.SpaceReducer.Eps = 1e-9;
            
            t1 = m.off1_createParamSamples;
            t2 = m.off2_genTrainingData;
            t3 = m.off3_computeReducedSpace;
            t4 = m.off4_genApproximationTrainData;
            
            factors = general.Utils.createCombinations([5 10 15 20], [1 2 3 4 5],[.1 .01 .001]);
            n = size(factors,2);
            K = approx.KernelApprox.empty;
            t5 = zeros(1,n);
            for i = 1:n
                fprintf('Starting approximation run %d of %d..\n',i,n);
                aa.NumGammas = factors(1,i);
                aa.MaxGFactor = factors(2,i);
                aa.MinGFactor = factors(3,i);
                t5(i) = m.off5_computeApproximation;
                K(i) = m.Approx.clone;
            end
            times = [t1 t2 t3 t4 t5];
            
            save nonzerox0_300s_dt01_largeAA_subsp m K times factors;
        end
        
        function MinMaxAdaptiveCWKA(T, dt)
            % Configuration to try and find the good old approximation with
            % little approximation error
            m = models.pcd.PCDModel(1);
            
            if nargin == 1
                dt = T/2000;
            elseif nargin == 0
                T = 6;
                dt = .001;
            end
            m.T = T; %[s]
            m.dt = dt; %[s]
            m.System.h = m.L/24;
            o = solvers.ode.MLode15i;
            o.AbsTol = 1e-6;
            o.RelTol = 1e-5;
            o.MaxStep = [];
            m.ODESolver = o;
            
            m.System.Params(1).Range = [1e-3 .1];
            m.System.Params(1).Desired = 12;
            
            s = sampling.GridSampler;
            s.Spacing = 'lin';
            m.Sampler = s;
            
            a = m.Approx;
%             s = approx.selection.DefaultSelector;
            %s = approx.selection.LinspaceSelector;
            s = approx.selection.TimeSelector;
            s.Size = 25000;
            a.TrainDataSelector = s;
            aa = approx.algorithms.MinMaxAdaptiveCWKA;
            aa.MaxExpansionSize = 80;
            aa.MaxRelErr = 1e-5;
            aa.MaxAbsErrFactor = 1e-5;
            aa.CheckMaxErrorPercent = .1;
            aa.InitialCenter = 't0';
            a.Algorithm = aa;
            
            m.Data = data.MemoryModelData;
            
            % Zero initial conditions
            %dim = size(m.System.x0.evaluate([]),1);
            %m.System.x0 = dscomponents.ConstInitialValue(zeros(dim,1));
            
            m.SpaceReducer = [];
%             m.SpaceReducer = spacereduction.PODGreedy;
%             m.SpaceReducer.Eps = 1e-9;
            
            a = KerMor.App;
            a.Verbose = 2;

            t1 = m.off1_createParamSamples;
            t2 = m.off2_genTrainingData;
            t3 = m.off3_computeReducedSpace;
            t4 = m.off4_genApproximationTrainData;
            
            %factors = general.Utils.createCombinations([5 10 15 20], [1 2 3 4 5],[.1 .01 .001]);
            factors = general.Utils.createCombinations(3, 10, 2);
            n = size(factors,2);
            K = approx.KernelApprox.empty;
            t5 = zeros(1,n);
            for i = 1:n
                fprintf('Starting approximation run %d of %d..\n',i,n);
                aa.NumGammas = factors(1,i);
                aa.MaxGFactor = factors(2,i);
                aa.MinGFactor = factors(3,i);
                t5(i) = m.off5_computeApproximation;
                K(i) = m.Approx.clone;
            end
            times = [t1 t2 t3 t4 t5];%#ok
            
            file = sprintf('MinMaxAdaptiveCWKA_T%d_dt%s', T, strrep(sprintf('%f',dt),'.','_'));
            save(file, 'm', 'K', 'times', 'factors');
        end
        
        function MinMaxAdaptiveCWKA_newrange(T, dt)
            % Configuration to try and find the good old approximation with
            % little approximation error
            m = models.pcd.PCDModel(1);
            
            if nargin == 1
                dt = T/2000;
            elseif nargin == 0
                T = 6;
                dt = .001;
            end
            m.T = T; %[s]
            m.dt = dt; %[s]
            m.System.h = m.L/24;
            o = solvers.ode.MLode15i;
            o.AbsTol = 1e-6;
            o.RelTol = 1e-5;
            o.MaxStep = [];
            m.ODESolver = o;
            
            m.System.Params(1).Range = [1e-5 .01];
            m.System.Params(1).Desired = 12;
            
            s = sampling.GridSampler;
            s.Spacing = 'lin';
            m.Sampler = s;
            
            a = m.Approx;
%             s = approx.selection.DefaultSelector;
            %s = approx.selection.LinspaceSelector;
            s = approx.selection.TimeSelector;
            s.Size = 25000;
            a.TrainDataSelector = s;
            aa = approx.algorithms.MinMaxAdaptiveCWKA;
            aa.MaxExpansionSize = 80;
            aa.MaxRelErr = 1e-5;
            aa.MaxAbsErrFactor = 1e-5;
            aa.CheckMaxErrorPercent = .1;
            aa.InitialCenter = 't0';
            a.Algorithm = aa;
            
            m.Data = data.MemoryModelData;
            
            % Zero initial conditions
            %dim = size(m.System.x0.evaluate([]),1);
            %m.System.x0 = dscomponents.ConstInitialValue(zeros(dim,1));
            
            m.SpaceReducer = [];
%             m.SpaceReducer = spacereduction.PODGreedy;
%             m.SpaceReducer.Eps = 1e-9;
            
            a = KerMor.App;
            a.Verbose = 2;

            t1 = m.off1_createParamSamples;
            t2 = m.off2_genTrainingData;
            t3 = m.off3_computeReducedSpace;
            t4 = m.off4_genApproximationTrainData;
            
            %factors = general.Utils.createCombinations([5 10 15 20], [1 2 3 4 5],[.1 .01 .001]);
            factors = general.Utils.createCombinations(3, 10, 2);
            n = size(factors,2);
            K = approx.KernelApprox.empty;
            t5 = zeros(1,n);
            for i = 1:n
                fprintf('Starting approximation run %d of %d..\n',i,n);
                aa.NumGammas = factors(1,i);
                aa.MaxGFactor = factors(2,i);
                aa.MinGFactor = factors(3,i);
                t5(i) = m.off5_computeApproximation;
                K(i) = m.Approx.clone;
            end
            times = [t1 t2 t3 t4 t5];%#ok
            
            file = sprintf('MinMaxAdaptiveCWKA_T%d_dt%s', T, strrep(sprintf('%f',dt),'.','_'));
            save(file, 'm', 'K', 'times', 'factors');
        end
        
        function MinMaxAdaptiveCWKA_newrange_des20(T, dt)
            % Configuration to try and find the good old approximation with
            % little approximation error
            m = models.pcd.PCDModel(1);
            
            if nargin == 1
                dt = T/2000;
            elseif nargin == 0
                T = 6;
                dt = .001;
            end
            m.T = T; %[s]
            m.dt = dt; %[s]
            m.System.h = m.L/24;
            o = solvers.ode.MLode15i;
            o.AbsTol = 1e-6;
            o.RelTol = 1e-5;
            o.MaxStep = [];
            m.ODESolver = o;
            
            m.System.Params(1).Range = [1e-5 .01];
            m.System.Params(1).Desired = 20;
            
            s = sampling.GridSampler;
            s.Spacing = 'lin';
            m.Sampler = s;
            
            a = m.Approx;
%             s = approx.selection.DefaultSelector;
            %s = approx.selection.LinspaceSelector;
            s = approx.selection.TimeSelector;
            s.Size = 25000;
            a.TrainDataSelector = s;
            aa = approx.algorithms.MinMaxAdaptiveCWKA;
            aa.MaxExpansionSize = 80;
            aa.MaxRelErr = 1e-5;
            aa.MaxAbsErrFactor = 1e-5;
            aa.CheckMaxErrorPercent = .1;
            aa.InitialCenter = 't0';
            a.Algorithm = aa;
            
            m.Data = data.MemoryModelData;
            
            % Zero initial conditions
            %dim = size(m.System.x0.evaluate([]),1);
            %m.System.x0 = dscomponents.ConstInitialValue(zeros(dim,1));
            
            m.SpaceReducer = [];
%             m.SpaceReducer = spacereduction.PODGreedy;
%             m.SpaceReducer.Eps = 1e-9;
            
            a = KerMor.App;
            a.Verbose = 2;

            t1 = m.off1_createParamSamples;
            t2 = m.off2_genTrainingData;
            t3 = m.off3_computeReducedSpace;
            t4 = m.off4_genApproximationTrainData;
            
            %factors = general.Utils.createCombinations([5 10 15 20], [1 2 3 4 5],[.1 .01 .001]);
            factors = general.Utils.createCombinations(3, 10, 2);
            n = size(factors,2);
            K = approx.KernelApprox.empty;
            t5 = zeros(1,n);
            for i = 1:n
                fprintf('Starting approximation run %d of %d..\n',i,n);
                aa.NumGammas = factors(1,i);
                aa.MaxGFactor = factors(2,i);
                aa.MinGFactor = factors(3,i);
                t5(i) = m.off5_computeApproximation;
                K(i) = m.Approx.clone;
            end
            times = [t1 t2 t3 t4 t5];%#ok
            
            file = sprintf('MinMaxAdaptiveCWKA_T%d_dt%s_newrange_des20', T, strrep(sprintf('%f',dt),'.','_'));
            save(file, 'm', 'K', 'times', 'factors');
        end
        
        function MinMaxAdaptiveCWKA_newR_des20log_subspa(T, dt)
            % Configuration to try and find the good old approximation with
            % little approximation error
            m = models.pcd.PCDModel(1);
            
            if nargin == 1
                dt = T/2000;
            elseif nargin == 0
                T = 6;
                dt = .001;
            end
            m.T = T; %[s]
            m.dt = dt; %[s]
            m.System.h = m.L/24;
            o = solvers.ode.MLode15i;
            o.AbsTol = 1e-6;
            o.RelTol = 1e-5;
            o.MaxStep = [];
            m.ODESolver = o;
            
            m.System.Params(1).Range = [1e-5 .01];
            m.System.Params(1).Desired = 20;
            
            s = sampling.GridSampler;
            s.Spacing = 'log';
            m.Sampler = s;
            
            a = m.Approx;
%             s = approx.selection.DefaultSelector;
            %s = approx.selection.LinspaceSelector;
            s = approx.selection.TimeSelector;
            s.Size = 25000;
            a.TrainDataSelector = s;
            aa = approx.algorithms.MinMaxAdaptiveCWKA;
            aa.MaxExpansionSize = 80;
            aa.MaxRelErr = 1e-5;
            aa.MaxAbsErrFactor = 1e-5;
            aa.CheckMaxErrorPercent = .07;
            aa.InitialCenter = 't0';
            a.Algorithm = aa;
            
            m.Data = data.MemoryModelData;
            
            % Zero initial conditions
            %dim = size(m.System.x0.evaluate([]),1);
            %m.System.x0 = dscomponents.ConstInitialValue(zeros(dim,1));
            
%            m.SpaceReducer = [];
            m.SpaceReducer = spacereduction.PODGreedy;
            m.SpaceReducer.Eps = 1e-9;
            
            a = KerMor.App;
            a.Verbose = 2;

            t1 = m.off1_createParamSamples;
            t2 = m.off2_genTrainingData;
            t3 = m.off3_computeReducedSpace;
            t4 = m.off4_genApproximationTrainData;
            
            %factors = general.Utils.createCombinations([5 10 15 20], [1 2 3 4 5],[.1 .01 .001]);
            factors = general.Utils.createCombinations(3, 10, 2);
            n = size(factors,2);
            K = approx.KernelApprox.empty;
            t5 = zeros(1,n);
            for i = 1:n
                fprintf('Starting approximation run %d of %d..\n',i,n);
                aa.NumGammas = factors(1,i);
                aa.MaxGFactor = factors(2,i);
                aa.MinGFactor = factors(3,i);
                t5(i) = m.off5_computeApproximation;
                K(i) = m.Approx.clone;
            end
            times = [t1 t2 t3 t4 t5];%#ok
            
            file = sprintf('MinMaxAdaptiveCWKA_T%d_dt%s_newR_des20log_subspa', T, strrep(sprintf('%f',dt),'.','_'));
            conf = object2str(m);
            save(file, 'm', 'K', 'times', 'factors', 'conf');
        end
        
        function MMACWKA_20logPar_SR(T, dt, dim)
            % Configuration to try and find the good old approximation with
            % little approximation error
            m = models.pcd.PCDModel(1);
            
            if nargin < 3
                dim = 100;
                if nargin == 1
                    dt = T/2000;
                elseif nargin == 0
                    T = 6;
                    dt = .001;
                end
            end
            m.T = T; %[s]
            m.dt = dt; %[s]
            m.System.h = m.L/(dim/4);
            o = solvers.ode.MLode15i;
            o.AbsTol = 1e-6;
            o.RelTol = 1e-5;
            o.MaxStep = [];
            m.ODESolver = o;
            
            m.System.Params(1).Range = [1e-6 .01];
            m.System.Params(1).Desired = 20;
            
            s = sampling.GridSampler;
            s.Spacing = 'log';
            m.Sampler = s;
            
            a = m.Approx;
%             s = approx.selection.DefaultSelector;
            %s = approx.selection.LinspaceSelector;
            s = approx.selection.TimeSelector;
            s.Size = 25000;
            a.TrainDataSelector = s;
            aa = approx.algorithms.MinMaxAdaptiveCWKA;
            aa.MaxExpansionSize = 80;
            aa.MaxRelErr = 1e-5;
            aa.MaxAbsErrFactor = 1e-5;
            aa.CheckMaxErrorPercent = .06;
            aa.InitialCenter = 't0';
            a.Algorithm = aa;
            
            m.Data = data.FileModelData(m);
            
            % Zero initial conditions
            %dim = size(m.System.x0.evaluate([]),1);
            %m.System.x0 = dscomponents.ConstInitialValue(zeros(dim,1));
            
%            m.SpaceReducer = [];
            m.SpaceReducer = spacereduction.PODGreedy;
            m.SpaceReducer.Eps = 1e-9;
            
            a = KerMor.App;
            a.Verbose = 2;

            t1 = m.off1_createParamSamples;
            t2 = m.off2_genTrainingData;
            t3 = m.off3_computeReducedSpace;
            t4 = m.off4_genApproximationTrainData;
            
            %factors = general.Utils.createCombinations([5 10 15 20], [1 2 3 4 5],[.1 .01 .001]);
            factors = general.Utils.createCombinations([4 12], 8, .1);
            n = size(factors,2);
            K = approx.KernelApprox.empty;
            t5 = zeros(1,n);
            for i = 1:n
                fprintf('Starting approximation run %d of %d..\n',i,n);
                aa.NumGammas = factors(1,i);
                aa.MaxGFactor = factors(2,i);
                aa.MinGFactor = factors(3,i);
                t5(i) = m.off5_computeApproximation;
                K(i) = m.Approx.clone;
            end
            times = [t1 t2 t3 t4 t5];%#ok
            
            file = sprintf('pcd1D_MMACWKA_20logPar_SR_T%d_dt%s_size%d', T, strrep(sprintf('%f',dt),'.','_'),dim);
            conf = object2str(m);
            save(file, 'm', 'K', 'times', 'factors', 'conf');
        end
    end
    
end