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
                    dim = 2000;
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
                this.System.MaxTimestep = this.dt;
                this.Name = '1D Burgers equation (A + f parts)';
            end
            this.Dimension = dim;
            
            this.SpaceReducer = spacereduction.PODGreedy;
            this.SpaceReducer.Eps = 1e-9;
            
            a = approx.DEIM;
            a.MaxOrder = 80;
            this.Approx = a;
            
            this.ErrorEstimator = error.DEIMEstimator;
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
        
        function m = test_Burgers_DEIM(dim, version)
            if nargin < 2
                version = 1;
                if nargin < 1
                    dim = 2000;
                end
            end
            m = models.burgers.Burgers(dim, version);
            
            if dim < 1000
                m.Data = data.MemoryModelData;
            else
                m.Data = data.FileModelData;
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
            eval(sprintf('save test_Burgers_DEIM_d%d_v%d',dim,version));
        end
        
        function m = test_Burgers_DEIM_B(dim)
            if nargin < 1
                dim = 200;
            end
            m = models.burgers.Burgers(dim, 2);
            
            if dim < 1000
                m.Data = data.MemoryModelData;
            else
                m.Data = data.FileModelData;
            end
            
            %% Sampling - manual
            s = sampling.ManualSampler;
            s.Samples = logspace(log10(0.04),log10(0.08),100);
            m.Sampler = s;
            
            %% Approx
            a = approx.DEIM;
            a.MaxOrder = 80;
            m.Approx = a;
            
            s = m.System;
            s.x0 = dscomponents.ConstInitialValue(zeros(dim,1));
            
            s.Inputs{1} = @(t)2*sin(2*t*pi);
            B = zeros(dim,1);
            B([4 5 6 7]) = [1 2 4 1];
            B(55:60) = 1;
            s.B = dscomponents.LinearInputConv(B);
            m.TrainingInputs = 1;
            
            offline_times = m.offlineGenerations;
            gitbranch = KerMor.getGitBranch;
            
            clear a s;
            eval(sprintf('save test_Burgers_DEIM_B_d%d',dim));
        end
    end
    
end