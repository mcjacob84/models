classdef WH10System < models.BaseDynSystem
    % Numerical experiments class for Paper WH10
    %
    % Current version works with KerMor 0.4
    
    properties(SetObservable)
        % SV's
        %
        % @propclass{experimental}
        svNum;
    end
    
    methods
        
        function this = WH10System(model)
            
            this = this@models.BaseDynSystem(model);
            
            this.registerProps('svNum');
            
            %% System settings
            this.MaxTimestep = [];
            
            % Sample bases
            this.svNum = 20;
            space = linspace(0,50,this.svNum);
            
            f = dscomponents.ParamTimeKernelCoreFun(this);
            fe = f.Expansion;
            fe.Centers.xi = repmat(space,model.dim,1);
            fe.Centers.ti = [];
            fe.Centers.mui = [];
            
            %% BaseCompWiseKernelApprox settings
            % Choose the kernel with such that each kernel vanishes (<
            % kerneleps) after kernelsupport support vectors.
            kernelsupport = 2;
            kerneleps = 0.00001;
            
            d = sqrt(model.dim)*space(2);
            %gamma = -((kernelsupport*d)^2)/log(kerneleps) preprint
            gamma = -((kernelsupport*d))/log(kerneleps)
            fe.Kernel = kernels.GaussKernel(gamma);
            fe.TimeKernel = kernels.NoKernel;
            fe.ParamKernel = kernels.NoKernel;
            this.f = f;
        end
    end
end

