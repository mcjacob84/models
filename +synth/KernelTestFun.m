classdef KernelTestFun < dscomponents.ParamTimeKernelCoreFun
    % KernelTestFun:
    %
    % Extracted these details from the former KernelTestSys.
    %
    % @author Daniel Wirtz @date 2011-07-07
    %
    % @new{0,5,dw,2011-07-07} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.morepas.org/software/index.html
    % - \c Documentation http://www.morepas.org/software/kermor/index.html
    % - \c License @ref licensing
    
    properties(SetObservable)
        % SV's
        %
        % @propclass{experimental} Test quantity.
        svNum;
        
        % @propclass{experimental} System dimension.
        dim;
    end
    
    methods
        function this = KernelTestFun(dim, pos_flag)
            
            this = this@dscomponents.ParamTimeKernelCoreFun;
            this.registerProps('svNum');
            
            if nargin == 1
                pos_flag = false;
            end
            
            % Sample bases
            this.svNum = 20;
            kexp = kernels.ParamTimeKernelExpansion;
            kexp.Centers.xi = repmat(linspace(-24,24,this.svNum),dim,1);
            kexp.Centers.ti = .5*(1:this.svNum);
            kexp.Centers.mui = rand(2,this.svNum);
            
            % Function coefficients
            offset = .5;
            %offset = 0;
            if pos_flag
                offset = 0;
            end
            ai = (rand(1,this.svNum)-offset);
            
            kexp.Ma = repmat(ai,dim,1);
            
            %% KernelApprox settings
            %this.Kernel = kernels.GaussKernel(4*sqrt(dims));
            kexp.Kernel = kernels.GaussKernel(sqrt(dim));
            kexp.Kernel.G = 1;
            kexp.TimeKernel = kernels.GaussKernel(.4);
            kexp.TimeKernel.G = 1;
            kexp.ParamKernel = kernels.GaussKernel(.4);
            kexp.ParamKernel.G = 1;
            this.Expansion = kexp;
            
            this.dim = dim;
        end
        
        function showBaseFun(this)
            % Debug method; displays the core function for each parameter
            % sample.
            figure;
            x = repmat(linspace(-50,50,max(this.svNum*5,100)),this.dim,1);
            fx = this.evaluate(x,0,[]);
            plot(x(1,:),fx(1,:),'r');
            xlabel('x'); ylabel('f(x)');
            title('KernelTest base function');
        end
    end
    
end