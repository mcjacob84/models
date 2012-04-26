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
        function this = Burgers(dim)
            if nargin < 1
                dim = 100;
            end
            this.T = 1;
            this.dt = .01;
            this.System = models.burgers.BurgersSys(this);
            this.Dimension = dim;
            this.ODESolver = solvers.ode.MLode15i;
            
            this.SpaceReducer = spacereduction.PODGreedy;
            this.SpaceReducer.Eps = 1e-8;
            
            a = this.Approx;
            a.ParamKernel = kernels.GaussKernel;
            al = approx.algorithms.VectorialKernelOMP;
            al.UseOGA = true;
            al.NumGammas = 30;
            al.MaxExpansionSize = 400;
            al.gameps = 1e-2;
            a.Algorithm = al;
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
            s.Samples = logspace(log10(0.005),log10(0.9),50);
            
            %% Approx
            a = approx.KernelApprox;
            a.Kernel = kernels.GaussKernel;
            a.ParamKernel = kernels.GaussKernel;
            
            al = approx.algorithms.VectorialKernelOMP;
            al.MaxGFactor = [1 0 1];
            al.MinGFactor = [.2 0 .6];
            al.gameps = 1e-4;
            al.MaxExpansionSize = 600;
            al.MaxAbsErrFactor = 1e-5;
            al.MaxRelErr = 1e-3;
            al.NumGammas = 25;
            a.Algorithm = al;
            m.Approx = a;
            
            m.offlineGenerations;
            save test_Burgers;
        end
    end
    
end