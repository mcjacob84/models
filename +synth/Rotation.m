classdef Rotation < models.BaseFullModel
    %ROTATION Synthetic rotation model
    
    methods
        function this = Rotation
            this.Verbose = 0;
            
            this.T = 25;
            this.dt = .1;
            
            %s = sampling.RandomSampler;
            %s.Samples = 10;
            %this.Sampler = s;
            this.Sampler = sampling.GridSampler;
            
            this.System = models.synth.RotationDynSys;
            this.Name = 'Synthetic rotation model';
            
            % Space reduction setup -> no reduction as only 2D!
%             sr = spacereduction.PODReducer;
%             sr.Mode = 'rel';
%             sr.Value = .2;
%             this.SpaceReducer = sr;
%             
            % Core Approximation
%             a = approx.CompWiseLS;
%             a.TimeKernel = kernels.GaussKernel;
%             a.SystemKernel = kernels.GaussKernel(2);
%             a.ParamKernel = kernels.LinearKernel;
%             a.lambda = 2;
            
%             a = approx.CompWiseSVR;
%             a.TimeKernel = kernels.GaussKernel;
%             a.SystemKernel = kernels.GaussKernel(2);
%             a.ParamKernel = kernels.LinearKernel;
%             a.eps = .05;
%             a.C = 100;
%             this.Approx = a;
            
            a = approx.CompWiseInt;
            a.TimeKernel = kernels.GaussKernel;
            a.SystemKernel = kernels.GaussKernel(2);
            %a.ParamKernel = kernels.LinearKernel;
            a.ParamKernel = kernels.GaussKernel(2);
            this.Approx = a;
            
            this.ODESolver = @ode45;
        end
        
    end
    
end
