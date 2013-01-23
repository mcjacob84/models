classdef RiemannBurgers < models.rbmatlab.RBMatlabModel
    %RIEMANNBURGERS Summary of this class goes here
    %   Detailed explanation goes here
            
    methods
        function this = RiemannBurgers(xnum, ynum, T)
            % Call superconstructor
            this = this@models.rbmatlab.RBMatlabModel;
            
            % Setup rbmatlab param struct
            params = struct;
            if nargin < 3
                T = .5;
                if nargin < 2
                    if nargin < 1
                         xnum = 100;    
                    end
                    ynum = 20;
                end
            end
            params.xnumintervals = xnum;
            params.ynumintervals = ynum;
            
            % Get rbmatlab model
            m = riemann_burgers_model(params);
            m.newton_regularisation = 0;
            m.T = T;
            m.nt = 100;
            % for a detailed simulation we do not need the affine parameter
            % decomposition
            m.decomp_mode = 0;
            this.RBMModel = m;
            
            this.RBMDataCont = models.rbmatlab.RBMDataContainer(gen_model_data(m));
            
            %this.ComputeParallel = false;
            
            % Times
            this.T = m.T;
            this.dt = m.T / m.nt;
            
            % System setup
            this.System = models.rbmatlab.RiemBurgSys_Fun(this);
            
            % Solver
            this.ODESolver = solvers.ExplEuler(this.dt);
            
            % Sampling
            this.Sampler = sampling.GridSampler;
            
            aa = approx.algorithms.FixedCompWiseKernelApprox;
            aa.ComputeParallel = false;
            a = approx.KernelApprox;
            a.Algorithm = aa;
            a.Kernel = kernels.GaussKernel(60);
            a.TimeKernel = kernels.GaussKernel(10*m.T/m.nt);
            a.ParamKernel = kernels.NoKernel;
            this.Approx = a;
            this.preApproximationTrainingCallback = @this.preApproxCallback;
            this.postApproximationTrainingCallback = @this.postApproxCallback;
            
            % Space reduction; choose only first row for subspace as the
            % problem repeats in y-direction.
            V = repmat(eye(params.xnumintervals),...
                params.ynumintervals,1)*sqrt(1/params.ynumintervals);
            this.SpaceReducer = spacereduction.ManualReduction(V,V);
        end
        
        function plot(this, t, y)%#ok
            sd.U = y;
            ppar.no_lines = 1;
            plot_sim_data(this.RBMModel, this.RBMDataCont.RBMData, sd, ppar);
        end
        
    end
    
    methods(Access=private)
        function preApproxCallback(this)
            xnum = this.RBMModel.xnumintervals;
            %this.Data.ApproxTrainData = this.Data.ApproxTrainData(1:xnum+3,:);
            this.Data.ApproxTrainData.fxi = this.Data.ApproxTrainData.fxi(1:xnum,:);
        end
        
        function postApproxCallback(this)
            ynum = this.RBMModel.ynumintervals;
            this.Data.ApproxTrainData.fxi = repmat(this.Data.ApproxTrainData.fxi,ynum,1);
            this.Approx.Ma = repmat(this.Approx.Ma,ynum,1);
            this.Approx.off = repmat(this.Approx.off,ynum,1);
        end
    end
    
end

