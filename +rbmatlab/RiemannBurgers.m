classdef RiemannBurgers < models.rbmatlab.BaseRBMatlabWrapper & models.BaseDynSystem & dscomponents.ACoreFun
    %RIEMANNBURGERS Summary of this class goes here
    %   Detailed explanation goes here
        
    methods
        function this = RiemannBurgers(xnum, ynum)
            % Call superconstructor
            this = this@models.rbmatlab.BaseRBMatlabWrapper;            
            
            % Setup rbmatlab param struct
            params.ynumintervals = 50;
            params.xnumintervals = 80;
            if nargin > 0
                params.xnumintervals = xnum;
                if nargin > 1
                    params.ynumintervals = ynum;        
                end
            end
            
            % Get rbmatlab model
            m = riemann_burgers_model(params);
            m.newton_regularisation = 0;
            
            d = models.rbmatlab.RBMDataContainer(gen_model_data(m));
            
            % Times
            %m.T = 0.004;
            %m.nt = 1;
            this.T = m.T;
            this.dt = m.T / m.nt;
            
            % Solver
            this.ODESolver = solvers.ExplEuler(this.dt);
            
            %this.PODFix.Value = 3;
            
            % Sampling
            this.Sampler = sampling.GridSampler;
            
            % Approximation: Use same kernels
            a = approx.CompWiseInt;
            k = kernels.GaussKernel(100);
            a.SystemKernel = k;
            a.TimeKernel = k;
            a.ParamKernel = k;
            this.Approx = a;
            
            % Space reduction; choose only first row for subspace as the
            % problem repeats in y-direction.
            V = repmat(eye(params.xnumintervals),...
                params.ynumintervals,1)*sqrt(1/params.xnumintervals);
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
            this.addParam('ULeft', [.1, .4], 3);
            this.addParam('URight', [.6, 1], 3);
            this.addParam('xFlux', [-1, -.1], 3);
            
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
            plot_sim_data(this.RBMModel, this.RBMDataCont.RBMData, sd, []);
        end
    end
    
    methods(Access=private)
        function x0 = getx0(this, mu) 
            this.RBMModel = this.RBMModel.set_mu(this.RBMModel,mu);
            % initial values by midpoint evaluation
            x0 = this.RBMModel.init_values_algorithm(this.RBMModel,this.RBMDataCont.RBMData);
        end
    end
    
end

