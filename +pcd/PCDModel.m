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
            
            this.T = 400; %[s]
            this.dt = .25; %[s]
            
            %s = sampling.RandomSampler;
            %s.Samples = 10;
            %this.Sampler = s;
            this.Sampler = sampling.GridSampler;
            
            switch dim
                case 2 
                    s = models.pcd.PCDSystem2D;
                    this.Name = 'Programmed Cell Death 2D';
                otherwise
                    s = models.pcd.PCDSystem1D;
                    this.Name = 'Programmed Cell Death 1D';
            end
            s.Model = this;
            s.updateSimConstants;
            this.System = s;
            
            % Space reduction setup
            sr = spacereduction.PODReducer;
            sr.Mode = 'rel';
            sr.Value = .3;
            %this.SpaceReducer = sr;
            this.SpaceReducer = [];
            
            % Core Approximation
%             a = approx.DefaultCompWiseKernelApprox;
%             a.CoeffComp = general.regression.KernelLS;
%             a.TimeKernel = kernels.GaussKernel;
%             a.SystemKernel = kernels.GaussKernel(2);
%             a.ParamKernel = kernels.LinearKernel;
%             a.lambda = 2;

            a = approx.AdaptiveCompWiseKernelApprox;
            a.MaxExpansionSize = 150;
            a.MaxRelErr = 1e-5;
            a.MaxAbsErrFactor = 1e-5;
            s = approx.selection.TimeSelector;
            s.MaxSize = 10000;
            a.TrainDataSelector = s;
            %a.ScalarSVR = general.regression.ScalarNuSVR;
            %a.ScalarSVR.nu = .6;
            %a.TimeKernel = kernels.LinearKernel;
            a.SystemKernel = kernels.GaussKernel;
            a.ParamKernel = kernels.GaussKernel;
            
%             a = approx.DefaultCompWiseKernelApprox;
%             a.CoeffComp = general.regression.ScalarEpsSVR;
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
    
%     methods (Static, Access=protected)
%         function obj = loadobj(s)
%             % Loads the properties for the PCDModel part of this
%             % class.
%             %
%             % See also: ALoadable BaseFullModel.loadobj
%             obj = models.pcd.PCDModel;
%             obj = loadobj@models.BaseFullModel(s, obj);
%         end
%     end
    
end

