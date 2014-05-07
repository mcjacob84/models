classdef KernelTestSys < models.BaseDynSystem
    % Kernel core function test model 1
    %
    % This class implements the dynamical system!
            
    methods
        function this = KernelTestSys(m, pos_flag)
            
            this = this@models.BaseDynSystem(m);
            
            if nargin < 2
                pos_flag = false;
            end
            
            this.addParam('mu1',[0 1],3);
            this.addParam('mu2',[0 .5],4);
            this.MaxTimestep = [];
            
            %% System settings
            this.x0 = dscomponents.ConstInitialValue(ones(this.Model.dim,1)*.5);
            
            this.f = models.synth.KernelTestFun(this, pos_flag);
        end
    end
end

