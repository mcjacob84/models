classdef PCDModel < models.BaseFullModel
% Base class for both 1D and 2D pcd models
%
% @author Daniel Wirtz @date 2011-03-16
%
% @change{0,6,dw,2011-11-27} Moved the static reduction experiments to their respective xD
% systems.
%
% @new{0,5,dw,2011-11-02} Added many static reduction experiment cases to the model.
%
% @new{0,3,dw,2011-03-16} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(Constant)
        % Typical diffusion rate for Caspase-8
        d1 = 1.8e-11; %[m^2/s]
        
        % Typical diffusion rate for Caspase-3
        d2 = 1.86e-11; %[m^2/s]
        
        % Typical diffusion rate for Pro-Caspase-8
        d3 = 1.89e-11; %[m^2/s]
        
        % Typical diffusion rate for Pro-Caspase-3
        d4 = 2.27e-11; %[m^2/s]
        
        % Typical cell length (from 1D)
        L = 1e-5; %[m]
    end
    
    methods
        function this = PCDModel(dim)
            % Creates a new instance of the PCDModel
            %
            % Parameters:
            % dim: The dimension to use @default 1
            
            % Use the PCDSystem
            if nargin == 0
                dim=1;
            end
            
            this.T = 6; %[s]
            this.dt = .001; %[s]
            % time scaling
            this.tau = this.L^2/this.d1;
            
            this.Data = data.FileModelData(this);
            
%             s = sampling.RandomSampler;
%             s.Samples = 10;
            s = sampling.GridSampler;
            s.Spacing = 'log';
            this.Sampler = s;
            
%             this.ODESolver = solvers.ode.MLWrapper(@ode23);
%             this.ODESolver = solvers.ode.ExplEuler(this.dt);
            o = solvers.ode.MLode15i;
            o.AbsTol = 1e-6;
            o.RelTol = 1e-5;
            o.MaxStep = [];
            this.ODESolver = o;
            
            switch dim
                case 3
                    s = models.pcd.PCDSystem3D(this);
                    this.Name = 'Programmed Cell Death 3D';
                case 2 
                    s = models.pcd.PCDSystem2D(this);
                    this.Name = 'Programmed Cell Death 2D';
                otherwise
                    s = models.pcd.PCDSystem1D(this);
                    this.Name = 'Programmed Cell Death 1D';
            end
            this.System = s;
            
            % Space reduction setup
            sr = spacereduction.PODGreedy;
            sr.Eps = 1e-10;
            this.SpaceReducer = sr;
%             this.SpaceReducer = [];
            
            % Core Approximation
%             a = approx.algorithms.DefaultCompWiseKernelApprox;
%             a.CoeffComp = general.regression.KernelLS;
%             a.TimeKernel = kernels.GaussKernel;
%             a.Kernel = kernels.GaussKernel(2);
%             a.ParamKernel = kernels.LinearKernel;
%             a.lambda = 2;

            a = approx.KernelApprox;
            a.TimeKernel = kernels.NoKernel;%kernels.GaussKernel(1);
            %a.TimeKernel.G = 1;
            a.Kernel = kernels.GaussKernel(1);
            a.Kernel.G = 1;
            a.ParamKernel = kernels.GaussKernel(1);
            a.ParamKernel.G = 1;
            
            s = approx.selection.LinspaceSelector;
            s.Size = 15000;
            a.TrainDataSelector = s;
            
            aa = approx.algorithms.VectorialKernelOMP;
            aa.MaxExpansionSize = 200;
            aa.MaxRelErr = 1e-5;
            aa.MaxAbsErrFactor = 1e-5;
            aa.MaxGFactor = 50;
            aa.MinGFactor = .01;
            aa.NumGammas = 40;
            a.Algorithm = aa; 
            this.Approx = a;
        end
        
        function plot(this, t, y)
            % Overrides standard method and forwards to the system's plot
            % function. (they are 1D and 2D)
            this.System.plot(this, t, y);
        end
    end
    
    methods(Static)
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
    end
end

