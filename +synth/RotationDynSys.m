classdef RotationDynSys < models.BaseDynSystem & dscomponents.ICoreFun
    %ROTATIONDYNSYS Synthetic 2D dynamical system with rotation
    %   Also implements the ICoreFun interface as the target function is
    %   quite simple.
        
    methods
        function this = RotationDynSys
            this.x0 = @(mu)[0;1];
            
            % Add params
            this.addParam('Rotation-speed', [pi/2-.5, pi/2+.5], 3);
            this.addParam('Angle offset', [-0.02,0.02], 3);
%             this.addParam('Rotation-speed', pi/2, 1);
%             this.addParam('Angle offset', 0, 1);
            
            % This class implements the ICoreFun interface!
            this.f = this;
        end
        
        function plot(this, model, t, y)
            plot(y(2,:),y(1,:),'r');
            plot@models.BaseDynSystem(this,model,t,y);
        end
        
        function fx = evaluate(this, x, t, mu)%#ok
            % Implements ICoreFun.evaluate
            a = mu(1);
            b = a+mu(2);
            %b = a+mu(2)+sin(t)/2;
            A = [cos(a) -sin(b); 
                 sin(a) cos(b)];
            fx = A*x;
        end
        
        function projected = project(this, V)
            % Implements ICoreFun(::IProjectable).project
            % no projection for this model as there are only two
            % dimensions.
            projected = this;
        end
    end
    
end

