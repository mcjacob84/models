classdef PCDModel < models.BaseFullModel
    % Base class for both 1D and 2D pcd models
    
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
            
            this.Verbose = 0;
            
            this.T = 500; %[s]
            this.dt = .2; %[s]
            
            %s = sampling.RandomSampler;
            %s.Samples = 10;
            %this.Sampler = s;
            this.Sampler = sampling.GridSampler;
            
            switch dim
                case 2 
                    this.System = models.pcd.PCDSystem2D(this);
                    this.Name = 'Programmed Cell Death 2D';
                otherwise
                    this.System = models.pcd.PCDSystem1D(this);
                    this.Name = 'Programmed Cell Death 1D';
            end
            
            % Space reduction setup
            sr = spacereduction.PODReducer;
            sr.Mode = 'rel';
            sr.Value = .3;
            this.SpaceReducer = sr;
            
            % Core Approximation
%             a = approx.CompWiseLS;
%             a.TimeKernel = kernels.GaussKernel;
%             a.SystemKernel = kernels.GaussKernel(2);
%             a.ParamKernel = kernels.LinearKernel;
%             a.lambda = 2;

            a = approx.CompWiseSVR;
            a.ApproxExpansionSize = 300;
            a.ScalarSVR = general.regression.ScalarNuSVR;
            a.ScalarSVR.nu = .6;
            a.TimeKernel = kernels.LinearKernel;
            a.SystemKernel = kernels.GaussKernel(50);
            a.ParamKernel = kernels.GaussKernel(50);
            
%             a = approx.CompWiseSVR;
%             a.TimeKernel = kernels.GaussKernel;
%             a.SystemKernel = kernels.GaussKernel(2);
%             a.ParamKernel = kernels.LinearKernel;
%             a.eps = .05;
%             a.C = 100;
            this.Approx = a;
            
            %this.ODESolver = solvers.MLWrapper(@ode23);
            this.ODESolver = solvers.ExplEuler(this.dt);
        end
        
        function plot(this, t, y)
            % Overrides standard method and forwards to the system's plot
            % function. (they are 1D and 2D)
            this.System.plot(this,t,y);
        end
    end
    
end

