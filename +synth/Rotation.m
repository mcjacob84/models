classdef Rotation < models.BaseFullModel
    
    methods
        function this = Rotation
            this.T = 10;
            this.dt = .025;
            
            %s = sampling.RandomSampler;
            %s.Samples = 10;
            %this.Sampler = s;
            this.Sampler = sampling.GridSampler;
            
            this.System = models.synth.RotationDynSys(this);
            this.Name = 'Synthetic rotation model';
            
            this.SpaceReducer = [];
%             
            % Core Approximation
%             a = approx.algorithms.Componentwise;
%             a.CoeffComp = general.regression.KernelLS;
%             a.TimeKernel = kernels.GaussKernel;
%             a.Kernel = kernels.GaussKernel(2);
%             a.ParamKernel = kernels.LinearKernel;
%             a.lambda = 2;
            
%             a = approx.algorithms.Componentwise;
%             a.CoeffComp = general.regression.ScalarEpsSVR;
%             a.TimeKernel = kernels.GaussKernel;
%             a.Kernel = kernels.GaussKernel(2);
%             a.ParamKernel = kernels.LinearKernel;
%             a.eps = .05;
%             a.C = 100;
%             this.Approx = a;
            
            a = approx.KernelApprox(this.System);
            a.Expansion.TimeKernel = kernels.GaussKernel;
            a.Expansion.Kernel = kernels.GaussKernel(2);
            %a.ParamKernel = kernels.LinearKernel;
            a.Expansion.ParamKernel = kernels.GaussKernel(2);
            this.Approx = a;
            
            this.ODESolver = solvers.MLWrapper(@ode45);
        end
        
        function plot(~, t, y)
            pm = PlotManager(false,1,2);
            ax = pm.nextPlot('3d','3D plot over time (z-axis)','x_1','x_2');
            plot3(ax,y(2,:),y(1,:),t,'r');
            ax = pm.nextPlot('1d','Plot over time','t','x_1,x_2');
            plot(ax,t,y(2,:),'r',t,y(1,:),'b');
            pm.LeaveOpen = true;
            pm.done;
        end
        
    end
    
    methods(Static)
        function test_Rotation_Simulation
            m = models.synth.Rotation;
            for k = 1:4
                mu = m.getRandomParam;
                [t,y] = m.simulate(mu);
                m.plot(t,y);
            end
        end
    end
    
end

