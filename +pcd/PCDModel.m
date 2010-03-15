classdef PCDModel < models.BaseFullModel
    %PCDMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function this = PCDModel
            this.Name = 'Programmed Cell Death';
            this.Verbose = 0;
            
            this.T = 1;
            this.dt = .1;
            
            %s = sampling.RandomSampler;
            %s.Samples = 10;
            %this.Sampler = s;
            this.Sampler = sampling.GridSampler;
            
            % Use the PCDSystem
            this.System = models.pcd.PCDSystem;
            
            % Space reduction setup
            sr = spacereduction.PODReducer;
            sr.Mode = 'eps';
            sr.Value = .1;
            this.SpaceReducer = sr;
            
            % Core Approximation
            a = approx.CompWiseSVR;
            a.eps = .3;
            a.C = 1;
            this.Approx = a;
            
            this.ODESolver = @ode45;
        end
    end
    
end

