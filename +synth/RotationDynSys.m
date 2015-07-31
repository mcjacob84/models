classdef RotationDynSys < models.BaseFirstOrderSystem
    %ROTATIONDYNSYS Synthetic 2D dynamical system with rotation
    %   Also implements the ACoreFun interface as the target function is
    %   quite simple.
        
    methods
        function this = RotationDynSys(model)
            this = this@models.BaseFirstOrderSystem(model);
            
            this.NumStateDofs = 2;
            this.updateDimensions;
            
            this.x0 = dscomponents.ConstInitialValue([0; 1]);
            
            % Add params
            this.addParam('Rotation-speed', pi/2, 'Range', [pi/2-.5, pi/2+.5], 'Desired', 3);
            this.addParam('Angle offset', 0, 'Range', [-0.02,0.02], 'Desired', 3);

            this.f = models.synth.RotationFun(this);
        end
    end
    
end

