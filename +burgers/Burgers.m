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
    
    properties
        % Azimuth and elevation to use for model plotting.
        PlotAzEl = [];
        SaveTag = 'models.burgers';
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
                this.Name = '1D unsteady Burgers equation (combined RHS)';
            elseif version == 2
                this.System = models.burgers.BurgersSys_A(this);
                this.System.MaxTimestep = this.dt;
                this.ODESolver = solvers.ode.SemiImplicitEuler(this);
                this.System.MaxTimestep = this.dt;
                this.Name = '1D unsteady Burgers equation';
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
            axis tight;
            title(sprintf('%s: dim=%d',this.Name,this.fDim));
            xlabel('x');
            ylabel('t');
            zlabel('y');
            rotate3d on;
            if ~isempty(this.PlotAzEl)
                view(this.PlotAzEl);
            end
        end
        
        function set.Dimension(this, value)
            this.fDim = value;
            this.System.newDim;
            this.Data.SimCache.clearTrajectories;
        end
        
        function dim = get.Dimension(this)
            dim = this.fDim;
        end
    end    
end