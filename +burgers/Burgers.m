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
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
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
    end
    
    methods
        function this = Burgers(dim, version)
            % Creates a new instance of the Burgers model.
            % 
            % Parameters:
            % dim: The dimension @type integer @default 2000
            % version: The version of the model to use. Possible choices
            % are 1 for a combined right hand side (diffusion &
            % nonlinearity in same function) or 2 for a linear part A
            % containing the diffusion and a nonlinear part f for the
            % quadratic term. @type integer @default 1
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
                this.ODESolver = solvers.MLode15i;
                this.SaveTag = sprintf('burgers_1D_d%d_combRHS',dim);
                this.Name = sprintf('1D-%dd unsteady Burgers equation (combined RHS)',dim);
            elseif version == 2
                this.System = models.burgers.BurgersSys_A(this);
                this.System.MaxTimestep = this.dt;
                this.ODESolver = solvers.SemiImplicitEuler(this);
                this.System.MaxTimestep = this.dt;
                this.SaveTag = sprintf('burgers_1D_d%d',dim);
                this.Name = sprintf('1D-%dd unsteady Burgers equation',dim);
            end
            this.Dimension = dim;
            
            this.SpaceReducer = spacereduction.PODGreedy;
            this.SpaceReducer.Eps = 1e-9;
            
            a = approx.DEIM(this.System);
            a.MaxOrder = 80;
            this.Approx = a;
            
            this.ErrorEstimator = error.DEIMEstimator;
        end
        
        function plot(this, t, y, pm_ax)
            if nargin < 4
                pm_ax = PlotManager;
                pm_ax.LeaveOpen = true;
            end
            nt = length(t);
            y  = [zeros(nt,1) y' zeros(nt,1)]; % add boundaries
            xx = linspace(this.Omega(1), this.Omega(2), this.fDim+2);

            if ~ishandle(pm_ax)
                pm_ax = pm_ax.nextPlot('burgers',sprintf('%s: dim=%d',this.Name,this.fDim),'x','t');
            end
            surfc(pm_ax,xx,t,y);
            shading interp;
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