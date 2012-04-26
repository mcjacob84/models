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
        
        function [f, ax] = plot(this, t, y)
            nt = length(t);
            y  = [zeros(nt,1) y' zeros(nt,1)]; % add boundaries
            xx = linspace(this.Omega(1), this.Omega(2), this.fDim+2);

            f = figure;
            ax = surfc(xx,t,y);
            shading interp;
            title(['Sol of Full System (FD):dim = ' num2str(this.fDim)]);
            xlabel('x');
            ylabel('t');
            zlabel('y');
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
            m.Approx.Kernel = kernels.PolyKernel(2);
            m.Approx.Algorithm.Dists = [zeros(2,10); linspace(.1,.9,10)];
            m.Approx.Algorithm.MaxExpansionSize = 50;
            m.offlineGenerations;
            
        end
    end
    
end