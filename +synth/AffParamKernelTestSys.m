classdef AffParamKernelTestSys < models.BaseDynSystem
    % Kernel core function test model using affine parametric initial values, input and output.
    %
    % This class implements the dynamical system!
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.morepas.org/software/index.html
    % - \c Documentation http://www.morepas.org/software/kermor/index.html
    % - \c License @ref licensing    
        
   
    
    methods
        
        function this = AffParamKernelTestSys(m, pos_flag)
            
            if nargin < 2
                pos_flag = false;
            end
            
            this = this@models.BaseDynSystem(m);
            
            %% System settings
            this.MaxTimestep = [];
            ai = dscomponents.AffineInitialValue;
            for i=1:4
                ai.addMatrix(['exp(' num2str(i) '*mu(1))+mu(2)'],ones(this.Model.dim,1)*i);
            end
            this.x0 = ai;
            
            this.addParam('mu1',[0 1],3);
            this.addParam('mu2',[0 .5],4);
            
            this.Inputs{1} = @(t)[.5*sin(t); 2*exp(-t)];
            b = dscomponents.AffLinInputConv;
            for i=1:4
                b.addMatrix(['.2*mu(1)*mu(2)/' num2str(i)],rand(this.Model.dim,2));
            end
            this.B = b;
            
            this.f = models.synth.KernelTestFun(this.Model.dim, pos_flag);
        end        
    end
end

