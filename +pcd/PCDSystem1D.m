classdef PCDSystem1D < models.pcd.BasePCDSystem
    %PCDSYSTEM1D Summary of this class goes here
    %   Detailed explanation goes here
    %
    % @author Daniel Wirtz @date 
    %
    % @change{0,2,dw,2011-03-09} Implemented the model according to the
    % Paper <a
    % href="http://www.simtech.uni-stuttgart.de/publikationen/prints.php?ID=285" 
    % target="_blank">Death wins against life in a spatially extended
    % apoptosis model</a>
    
    methods
        function this = PCDSystem1D(model)
            this = this@models.pcd.BasePCDSystem(model);

            % Set core function
            this.f = models.pcd.CoreFun1D(this);
            
            % Spatial resolution (in real sizes)
            this.Omega = [0 1] * this.Model.L;
            this.h = this.Model.L/24;
           
            % Add input param (is getting inserted after the BasePCDSystem
            % condtructor params, so number 9!)
            this.addParam('U', [1e-5, 1e-2], 50);
        end

        function plot(this, model, t, y)
            % Performs a plot for this model's results.
            %
            
            if length(t) > 700
                idx = round(linspace(1,length(t),700));
                t = t(idx);
                y = y(:,idx);
            end
            states = {'alive','unstable','dead'};
            
            figure;
            X = t;
            Y = this.Omega(1):this.h:this.Omega(2);
            m = model.System.Dims(1);
            doplot(y(1:m,:),'Caspase-8 (x_a)',1);
            doplot(y(m+1:2*m,:),'Caspase-3 (y_a)',2);
            doplot(y(2*m+1:3*m,:),'Pro-Caspase-8 (x_i)',3);
            doplot(y(3*m+1:end,:),'Pro-Caspase-3 (y_i)',4);
            
            function doplot(y,thetitle,pnr)
                subplot(2,2,pnr);
                mesh(X,Y,y);
                %surf(X,Y,y,'EdgeColor','none');
                xlabel('Time [s]');
                ylabel(sprintf('%.2e to %.2e: Cell core to hull [m]',this.Omega(1),this.Omega(2)));
                grid off;
                mi = min(y(:));
                Ma = max(y(:));
                if abs((mi-Ma) / mi) < 1e-14
                    mi = .999*mi; Ma=1.001*Ma;
                end
                axis([0 this.Model.T this.Omega mi Ma]);
                di = abs(this.SteadyStates(:,pnr)-y(end));
                reldi = di ./ (this.SteadyStates(:,pnr)+eps);
                reldistr = general.Utils.implode(reldi,', ','%2.3e');
                if any(reldi > .1) || any(reldi < 10)
                    [~, id] = min(di);
                    title(sprintf('Model "%s", %s concentrations\nCell state at T=%d: %s\n%s', model.Name, thetitle,...
                    model.T,states{id},reldistr));
                else
                    title(sprintf('Model "%s", %s concentrations\n%s', model.Name, thetitle,reldistr));
                end
            end
        end
    end
    
    methods(Access=protected)
        function newSysDimension(this)
            % Assign fitting initial value
            m = prod(this.Dims);
            
            x0 = zeros(4*m,1);
            x0(1:2*m) = 1e-16;
            x0(2*m+1:end) = 1e-9;
            this.x0 = dscomponents.ConstInitialValue(x0);
            
            e = ones(m,1);
            A = spdiags([e -2*e e],-1:1,m,m)/this.hs^2;
            A(1,2) = 2/this.hs^2;
            A(end,end-1) = 2/this.hs^2;
            A = blkdiag(A,this.Diff(1)*A,this.Diff(2)*A,this.Diff(3)*A);
            this.A = dscomponents.LinearCoreFun(A);
            
            % Extracts the caspase-3 concentrations from the result
%             C = zeros(m,4*m);
%             C(:,m+1:2*m) = diag(ones(m,1));
            %this.C = dscomponents.LinearInputConv(sparse(C));
        end
    end
    
    methods(Static)
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

