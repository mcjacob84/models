classdef KernelTest < models.BaseFullModel
    % Kernel core function test model 1
    %
    % This class implements both the model and dynamical system!
    %
    
    properties
        % The system's dimension
        dim;
    end
    
    methods
        
        function this = KernelTest(dims, pos_flag)
            
            if nargin < 2
                pos_flag = false;
                if nargin < 1
                    dims = 1000;
                end
            end
            this.dim = dims;
            
            
            %% Model settings
            this.Name = 'Kernel test model';
            
            this.T = 5;
            this.dt = .08;
            
            this.Sampler = [];%sampling.GridSampler;
            
            % This class implements a fake Approx subclass to allow access
            % to the this.Ma property for the error estimator.
            this.Approx = [];
            
            s = spacereduction.PODReducer;
            s.UseSVDS = true;
            s.Mode = 'abs';
            s.Value = 1;
            this.SpaceReducer = s;
            
            %this.ODESolver = solvers.ode.MLWrapper(@ode45);
            this.ODESolver = solvers.ode.ExplEuler;
            %this.ODESolver = solvers.ode.Heun;
            
            %% System settings
            this.System = models.synth.KernelTestSys(this, pos_flag);
        end
        
    end
      
    methods(Static)
        
        function r = runTest(model)
            model.offlineGenerations;
            r = model.buildReducedModel;
            r.analyze;
        end
        
        function m = getTest1(varargin)
            m = models.synth.KernelTest(varargin{:});
            
            V = ones(m.dim,1)*sqrt(1/m.dim);
            m.SpaceReducer = spacereduction.ManualReduction(V);
        end
        
        function m = getTest2(varargin)
            m = models.synth.KernelTest(varargin{:});
            
            m.System.Inputs{1} = @(t)4;
            m.System.B = dscomponents.LinearInputConv(ones(m.dim,1));
            
            V = ones(m.dim,1)*sqrt(1/m.dim);
            m.SpaceReducer = spacereduction.ManualReduction(V);
        end
        
        function m = getTest3(varargin)
            m = models.synth.KernelTest(varargin{:});
            
            x0 = rand(m.dim,1);
            m.System.x0 = @(mu)x0;
        end
        
        function m = getTest4(varargin)
            m = models.synth.KernelTest(varargin{:});
            
            x0 = rand(m.dim,1);
            m.System.x0 = @(mu)x0;
            V = ones(m.dim,1)*sqrt(1/m.dim);
            m.SpaceReducer = spacereduction.ManualReduction(V);
        end
        
        function m = getTest5(varargin)
            m = models.synth.KernelTest(varargin{:});
            
            m.System.B = dscomponents.LinearInputConv(rand(m.dim,1));
            m.System.Inputs{1} = @(t)4;
        end
        
        function m = getTest6(varargin)
            m = models.synth.KernelTest(varargin{:});
            
            m.System.B = dscomponents.LinearInputConv(rand(m.dim,1));
            m.System.Inputs{1} = @(t)4;
            
            V = ones(m.dim,1)*sqrt(1/m.dim);
            m.SpaceReducer = spacereduction.ManualReduction(V);
        end
        
        function m = getTest7(varargin)
            m = models.synth.KernelTest(varargin{:});
            
            m.System.Inputs{1} = @(t)4;
            
            B = ones(m.dim,1);
            B(1:m.dim/2) = -1;
            m.System.B = dscomponents.LinearInputConv(B);
        end
        
        function m = getTest8(varargin)
            m = models.synth.KernelTest(varargin{:});
            
            m.System.Inputs{1} = @(t)4;
            
            B = ones(m.dim,1);
            B(1:m.dim/2) = -1;
            m.System.B = dscomponents.LinearInputConv(B);
            
            V = ones(m.dim,1)*sqrt(1/m.dim);
            m.SpaceReducer = spacereduction.ManualReduction(V);
        end
        
        function m = getTest9(varargin)
            m = models.synth.KernelTest(varargin{:});
            
            m.System.Inputs{1} = @(t)4;
            
            x0 = (rand(m.dim,1)-.5)*3;
            m.System.x0 = @(mu)x0;
            
            B = ones(m.dim,1);
            B(1:m.dim/2) = -1;
            m.System.B = dscomponents.LinearInputConv(B);
            
            V = ones(m.dim,1)*sqrt(1/m.dim);
            m.SpaceReducer = spacereduction.ManualReduction(V);
        end
        
        function m = getTest10(varargin)
            m = models.synth.KernelTest(varargin{:});
            m.T = 20;
            
            m.System.Inputs{1} = @(t)sin(2*t);
            
            x0 = (rand(m.dim,1)-.5)*3;
            m.System.x0 = @(mu)x0;
            
            B = ones(m.dim,1);
            B(1:m.dim/2) = -1;
            m.System.B = dscomponents.LinearInputConv(B);
            
            V = ones(m.dim,1)*sqrt(1/m.dim);
            m.SpaceReducer = spacereduction.ManualReduction(V);
        end
        
        function m = getTest11(varargin)
            m = models.synth.KernelTest(varargin{:});
            m.T = 20;
            
            m.System.B = dscomponents.LinearInputConv(rand(m.dim,1));
            m.System.Inputs{1} = @(t)sin(2*t);
        end
    end
    
end

