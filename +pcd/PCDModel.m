classdef PCDModel < models.BaseFullModel
    %PCDMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function this = PCDModel(dim)
            % Creates a new instance of the PCDModel
            this.Verbose = 0;
            
            this.T = 200;
            this.dt = .2;
            
            %s = sampling.RandomSampler;
            %s.Samples = 10;
            %this.Sampler = s;
            this.Sampler = sampling.GridSampler;
            
            % Use the PCDSystem
            if nargin == 0
                dim=1;
            end
            switch dim
                case 2 
                    this.System = models.pcd.PCDSystem2D;
                    this.Name = 'Programmed Cell Death 2D';
                otherwise
                    this.System = models.pcd.PCDSystem1D;
                    this.Name = 'Programmed Cell Death 1D';
            end
            
            % Space reduction setup
            sr = spacereduction.PODReducer;
            sr.Mode = 'rel';
            sr.Value = .2;
            this.SpaceReducer = sr;
            
            % Core Approximation
%             a = approx.CompWiseLS;
%             a.TimeKernel = kernels.RBFKernel;
%             a.SystemKernel = kernels.RBFKernel(2);
%             a.ParamKernel = kernels.LinearKernel;
%             a.lambda = 2;
            a = approx.CompWiseInt;
            a.TimeKernel = kernels.RBFKernel(2);
            a.SystemKernel = kernels.RBFKernel(2);
            a.ParamKernel = kernels.RBFKernel(2);
            
%             a = approx.CompWiseSVR;
%             a.TimeKernel = kernels.RBFKernel;
%             a.SystemKernel = kernels.RBFKernel(2);
%             a.ParamKernel = kernels.LinearKernel;
%             a.eps = .05;
%             a.C = 100;
            this.Approx = a;
            
            this.ODESolver = @ode23;
        end
    end
    
end

