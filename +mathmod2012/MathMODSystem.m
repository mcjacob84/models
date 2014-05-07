classdef MathMODSystem < models.BaseDynSystem
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
        
        function this = MathMODSystem(model)
            
            this = this@models.BaseDynSystem(model);
            
            this.registerProps('svNum');
            
            %% System settings
            dims = model.dim;
            
            this.MaxTimestep = [];
            
            this.addParam('input shift',[1 1],1);
            this.addParam('expansion param',[0 10],1);
            this.addParam('initial value',[0 1],1);
            this.DependentParamIndices = 2;
            
            % Sample bases
            this.svNum = 20;

            f = dscomponents.ParamTimeKernelCoreFun(this);
            fe = f.Expansion;
            this.f = f;
            
            % Choose the kernel with such that each kernel vanishes (<
            % kerneleps) after kernelsupport support vectors. 
            kernelsupport = 2;
            kerneleps = 0.00001;
            
            % Space kernel
            space = linspace(0,50,this.svNum);
            fe.Centers.xi = repmat(space,dims,1);
            d = sqrt(dims)*space(2);
            fe.Kernel = kernels.GaussKernel;
            fe.Kernel.setGammaForDistance(kernelsupport*d,kerneleps)
            fe.Kernel.G = 1;
            
            % Time kernel
            fe.TimeKernel = kernels.NoKernel;
            fe.Centers.ti = [];
            
            % Param kernel
            pk = kernels.GaussKernel;
            pspace = linspace(0,10,this.svNum);
            % Only the 2nd entry is used, so fill rest with zeros
            pk.P = 2;
            pk.G = 1;
            pk.setGammaForDistance(sqrt(this.ParamCount)*pspace(2)*20,kerneleps);
            fe.ParamKernel = pk;
            fe.Centers.mui = [zeros(1,this.svNum); pspace; zeros(1,this.svNum)];
        end
    end    
end

