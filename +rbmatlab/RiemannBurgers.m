classdef RiemannBurgers < models.rbmatlab.BaseRBMatlabWrapper & models.BaseDynSystem & dscomponents.ACoreFun
    %RIEMANNBURGERS Summary of this class goes here
    %   Detailed explanation goes here
            
    methods
        function this = RiemannBurgers(xnum, ynum, T)
            % Call superconstructor
            this = this@models.rbmatlab.BaseRBMatlabWrapper;            
            
            % Setup rbmatlab param struct
            params.ynumintervals = 20;
            params.xnumintervals = 100;
            if nargin > 0
                params.xnumintervals = xnum;
                if nargin > 1
                    params.ynumintervals = ynum;
                    if nargin == 2
                        T = 0.5;
                    end
                end
            end
            
            % Get rbmatlab model
            m = riemann_burgers_model(params);
            m.newton_regularisation = 0;
            m.T = T;
            m.nt = 100;
            
            d = models.rbmatlab.RBMDataContainer(gen_model_data(m));
            
            %this.ComputeParallel = false;
            
            % Times
            this.T = m.T;
            this.dt = m.T / m.nt;
            
            % Solver
            this.ODESolver = solvers.ExplEuler(this.dt);
            
            % Sampling
            this.Sampler = sampling.GridSampler;
            
            % Approximation: Use same kernels
            if false
                a = approx.DefaultCompWiseKernelApprox;
                svr = general.regression.ScalarNuSVR;
                svr.nu = .6;
                qp = solvers.qpOASES;
                svr.QPSolver = qp;
                a.CoeffComp = svr;
            else
                a = approx.DefaultCompWiseKernelApprox;
            end
            a.SystemKernel = kernels.GaussKernel(60);
            a.TimeKernel = kernels.GaussKernel(10*m.T/m.nt);
            a.ParamKernel = kernels.NoKernel;
            a.ComputeParallel = false;
            this.Approx = a;
            this.ApproxExpansionSize = 300;
            this.preApproximationTrainingCallback = @this.preApproxCallback;
            this.postApproximationTrainingCallback = @this.postApproxCallback;
            
            % Space reduction; choose only first row for subspace as the
            % problem repeats in y-direction.
            V = repmat(eye(params.xnumintervals),...
                params.ynumintervals,1)*sqrt(1/params.ynumintervals);
            this.SpaceReducer = spacereduction.ManualReduction(V,V);
            
            %% System setup
            this.System = this;
            
            % for a detailed simulation we do not need the affine parameter
            % decomposition
            m.decomp_mode = 0;
            this.x0 = @(mu)this.getx0(mu);
            
            % DS-Components
            this.f = this;
            this.B = [];
            this.Inputs = {};
            
            % Parameters
            this.addParam('ULeft', [.3, .3], 1);
            this.addParam('URight', [.1, 1], 30);
            this.addParam('xFlux', [-.5, -.5], 1);
            
            this.RBMModel = m;
            this.RBMDataCont = d;
            
            % Output conversion
            %c = zeros(1, params.xnumintervals*params.ynumintervals);
            %c(round(.75*params.xnumintervals)) = 1;
            %this.C = dscomponents.PointerOutputConv(@(t,mu)c,false);
        end
        
        function y = evaluateCoreFun(this, x, t, mu)
            this.RBMModel = this.RBMModel.set_mu(this.RBMModel,mu);
            this.RBMModel = this.RBMModel.set_time(this.RBMModel,t);
            y = -this.RBMModel.L_E_local_ptr(this.RBMModel, this.RBMDataCont.RBMData, x, []);
        end        
        
        function plot(this, t, y)
            sd.U = y;
            ppar.no_lines = 1;
            plot_sim_data(this.RBMModel, this.RBMDataCont.RBMData, sd, ppar);
        end
        
       
    end
    
    methods(Access=private)
        function x0 = getx0(this, mu)
            this.RBMModel = this.RBMModel.set_mu(this.RBMModel,mu);
            % initial values by midpoint evaluation
            x0 = this.RBMModel.init_values_algorithm(this.RBMModel,this.RBMDataCont.RBMData);
        end
        
        function preApproxCallback(this)
            xnum = this.RBMModel.xnumintervals;
            %this.Data.ApproxTrainData = this.Data.ApproxTrainData(1:xnum+3,:);
            this.Data.ApproxfValues = this.Data.ApproxfValues(1:xnum,:);
        end
        
        function postApproxCallback(this)
            ynum = this.RBMModel.ynumintervals;
            this.Data.ApproxfValues = repmat(this.Data.ApproxfValues,ynum,1);
            this.Approx.Ma = repmat(this.Approx.Ma,ynum,1);
            this.Approx.off = repmat(this.Approx.off,ynum,1);
        end
    end
    
end

