classdef Burgers < models.BaseFullModel
% Burgers: 
%
%
%
% @author Daniel Wirtz @date 2012-04-24
%
% @new{0,6,dw,2012-04-24} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(Constant)
        Omega = [0 1];
    end

    properties(Dependent)
        Dimension;
    end
    
    properties(Access=private)
        fDim;
    end
    
    methods
        function this = Burgers(dim, version)
            this = this@models.BaseFullModel;
            if nargin < 2
                version = 1;
                if nargin < 1
                    dim = 100;
                end
            end
            this.T = 1;
            this.dt = .01;
            if version == 1
                this.System = models.burgers.BurgersSys(this);
                this.ODESolver = solvers.ode.MLode15i;
                this.Name = '1D Burgers equation (combined RHS)';
            elseif version == 2
                this.System = models.burgers.BurgersSys_A(this);
                this.System.MaxTimestep = this.dt;
                this.ODESolver = solvers.ode.SemiImplicitEuler(this);
                this.Name = '1D Burgers equation (A + f parts)';
            end
            this.Dimension = dim;
            
            this.SpaceReducer = spacereduction.PODGreedy;
            this.SpaceReducer.Eps = 1e-8;
            
            a = approx.DEIM;
            a.MaxOrder = 50;
            this.Approx = a;
%             a.ParamKernel = kernels.GaussKernel;
%             al = approx.algorithms.VectorialKernelOMP;
%             al.UseOGA = true;
%             al.NumGammas = 30;
%             al.MaxExpansionSize = 400;
%             al.gameps = 1e-2;
%             a.Algorithm = al;
        end
        
        function [f, ax] = plot(this, t, y, ax)
            if nargin < 4
                f = figure;
                ax = gca(f);
            end
            nt = length(t);
            y  = [zeros(nt,1) y' zeros(nt,1)]; % add boundaries
            xx = linspace(this.Omega(1), this.Omega(2), this.fDim+2);

            surfc(ax,xx,t,y);
            shading interp;
            title(['Sol of Full System (FD):dim = ' num2str(this.fDim)]);
            xlabel('x');
            ylabel('t');
            zlabel('y');
            rotate3d on;
        end
        
        function set.Dimension(this, value)
            this.fDim = value;
            this.System.newDim;
            this.simCache.clearTrajectories;
        end
        
        function dim = get.Dimension(this)
            dim = this.fDim;
        end
    end
    
    methods(Static)
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
        
        function m = test_Burgers_DEIM
            m = models.burgers.Burgers;
            
            %% Sampling - manual
            s = sampling.ManualSampler;
            s.Samples = logspace(log10(0.04),log10(0.08),100);
            m.Sampler = s;
            
            %% Approx
            a = approx.DEIM;
            a.MaxOrder = 60;
            m.Approx = a;
            
%             m.off1_createParamSamples;
%             m.off2_genTrainingData;
%             m.off3_computeReducedSpace;
%             m.off4_genApproximationTrainData;
%             m.off5_computeApproximation;
            m.offlineGenerations;
            save test_Burgers_DEIM;
        end
    end
    
end