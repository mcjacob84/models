classdef KernelTestSys < models.BaseFirstOrderSystem
    % Kernel core function test model 1
    %
    % This class implements the dynamical system!
            
    methods
        function this = KernelTestSys(m, pos_flag)
            
            this = this@models.BaseFirstOrderSystem(m);
            
            this.NumStateDofs = this.Model.dim;
            this.updateDimensions;
            
            if nargin < 2
                pos_flag = false;
            end
            
            this.addParam('mu1', .5 , 'Range', [0 1], 'Desired', 3);
            this.addParam('mu2', .25, 'Range', [0 .5], 'Desired', 4);
            this.MaxTimestep = m.dt;
            
            %% System settings
            this.x0 = dscomponents.ConstInitialValue(ones(this.Model.dim,1)*.5);
            
            this.f = models.synth.KernelTestFun(this, pos_flag);
        end
    end
end

