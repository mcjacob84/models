classdef KernelTestSys < models.BaseDynSystem & dscomponents.CompwiseKernelCoreFun
    % Kernel core function test model 1
    %
    % This class implements the dynamical system!
        
    properties(SetObservable)
        % SV's
        %
        % @propclass{experimental} Test quantity.
        svNum;
    end
    
    methods
        
        function this = KernelTestSys(m, pos_flag)
            
            this = this@models.BaseDynSystem(m);
            this.registerProps('svNum');
            
            if nargin < 2
                pos_flag = false;
            end
            
            %% System settings
            this.x0 = @(mu)ones(this.Model.dim,1)*.5;
            
            this.f = this;
            
            %this.Inputs{1} = @(t)0;
            
            % Sample bases
            this.svNum = 20;
            this.Centers.xi = repmat(linspace(-24,24,this.svNum),this.Model.dim,1);
            %this.Centers.xi = linspace(-20,20,this.svNum);
            this.Centers.ti = [];
            this.Centers.mui = [];
            
            % Function coefficients
            offset = .5;
            %offset = 0;
            if pos_flag
                offset = 0;
            end
            ai = (rand(1,this.svNum)-offset);
            
            this.Ma = repmat(ai,this.Model.dim,1);
            
            %% BaseCompWiseKernelApprox settings
            %this.SystemKernel = kernels.GaussKernel(4*sqrt(dims));
            this.SystemKernel = kernels.GaussKernel(sqrt(this.Model.dim));
            this.TimeKernel = kernels.NoKernel;
            this.ParamKernel = kernels.NoKernel;
        end
        
        function c = getcfi(this, z, C, t, mu)
            C = 100;
            di = this.xi - repmat(this.Data.V*z,1,this.sv);
            di = sqrt(sum(di.^2));
            case1 = di - C >= 0;
            case2 = ~case1;
            b = this.sk.Gamma;
            t2 = exp(-((di+C).^2/b));
            ci(case1) = (exp(-(di(case1)-C).^2/b) - t2(case1));
            ci(case2) = (1 - t2(case2));
            c = this.Ma_norms*ci';
        end
        
        function showBaseFun(this)
            % Debug method; displays the core function for each parameter
            % sample.
            figure;
            x = repmat(linspace(-50,50,max(this.svNum*5,100)),this.Model.dim,1);
            fx = this.evaluate(x,0,[]);
            plot(x(1,:),fx(1,:),'r');
            xlabel('x'); ylabel('f(x)');
            title('KernelTest base function');
        end
        
    end
end

