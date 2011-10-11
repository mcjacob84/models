classdef PCDModel < models.BaseFullModel
    % Base class for both 1D and 2D pcd models
    
    properties(Constant)
        % Typical diffusion rate for Caspase-8
        d1 = 1.8e-11; %[m^2/s]
        
        % Typical diffusion rate for Caspase-3
        d2 = 1.86e-11; %[m^2/s]
        
        % Typical diffusion rate for Pro-Caspase-8
        d3 = 1.89e-11; %[m^2/s]
        
        % Typical diffusion rate for Pro-Caspase-3
        d4 = 2.27e-11; %[m^2/s]
        
        % Typical cell length (from 1D)
        L = 1e-5; %[m]
    end
    
    methods
        function this = PCDModel(dim)
            % Creates a new instance of the PCDModel
            %
            % Parameters:
            % dim: The dimension to use @default 1
            
            % Use the PCDSystem
            if nargin == 0
                dim=1;
            end
            
            this.T = 1; %[s]
            this.dt = .001; %[s]
            % time scaling
            this.tau = this.L^2/this.d1;
            
            this.Data = data.FileModelData(this);
            
            %s = sampling.RandomSampler;
            %s.Samples = 10;
            %this.Sampler = s;
            this.Sampler = sampling.GridSampler;
            
%             this.ODESolver = solvers.ode.MLWrapper(@ode23);
            this.ODESolver = solvers.ode.ExplEuler(this.dt);
%             o = solvers.ode.MLode15i;
%             o.AbsTol = 1e-5;
%             o.RelTol = 1e-5;
%             o.MaxStep = this.dt;
%             this.ODESolver = o;
            
            switch dim
                case 3
                    s = models.pcd.PCDSystem3D(this);
                    this.Name = 'Programmed Cell Death 3D';
                case 2 
                    s = models.pcd.PCDSystem2D(this);
                    this.Name = 'Programmed Cell Death 2D';
                otherwise
                    s = models.pcd.PCDSystem1D(this);
                    this.Name = 'Programmed Cell Death 1D';
            end
            s.MaxTimestep = this.dt;
            this.System = s;
            
            % Space reduction setup
            sr = spacereduction.TrajectoryGreedy;
            sr.Eps = 1e-6;
            this.SpaceReducer = sr;
            
            % Core Approximation
%             a = approx.algorithms.DefaultCompWiseKernelApprox;
%             a.CoeffComp = general.regression.KernelLS;
%             a.TimeKernel = kernels.GaussKernel;
%             a.Kernel = kernels.GaussKernel(2);
%             a.ParamKernel = kernels.LinearKernel;
%             a.lambda = 2;

            a = approx.KernelApprox;
            a.TimeKernel = kernels.GaussKernel(1);
            a.Kernel = kernels.GaussKernel(1);
            a.ParamKernel = kernels.GaussKernel(1);
            
            s = approx.selection.LinspaceSelector;
            s.Size = 10000;
            a.TrainDataSelector = s;
            
            aa = approx.algorithms.AdaptiveCompWiseKernelApprox;
            aa.MaxExpansionSize = 164;
            aa.MaxRelErr = 1e-5;
            aa.MaxAbsErrFactor = 1e-5;
            a.Algorithm = aa; 

            this.Approx = a;
            
            %a.ScalarSVR = general.regression.ScalarNuSVR;
            %a.ScalarSVR.nu = .6;
            
        end
        
        function plot(this, t, y)
            % Overrides standard method and forwards to the system's plot
            % function. (they are 1D and 2D)
            this.System.plot(this,t,y);
        end
    end
    
%     methods (Static, Access=protected)
%         function obj = loadobj(obj)
%             % Loads the properties for the PCDModel part of this
%             % class.
%             %
%             % See also: ALoadable BaseFullModel.loadobj
%             %
%             %obj = models.pcd.PCDModel;
%             %obj = loadobj@models.BaseFullModel(s, obj);
%         end
%     end
    
end

