classdef RotationDynSys < models.BaseDynSystem & dscomponents.ACoreFun
    %ROTATIONDYNSYS Synthetic 2D dynamical system with rotation
    %   Also implements the ACoreFun interface as the target function is
    %   quite simple.
        
    methods
        function this = RotationDynSys
            this = this@dscomponents.ACoreFun;
            this.TimeDependent = false;
            
            this.x0 = @(mu)[0;1];
            
            % Add params
            this.addParam('Rotation-speed', [pi/2-.5, pi/2+.5], 3);
            this.addParam('Angle offset', [-0.02,0.02], 3);
%             this.addParam('Rotation-speed', pi/2, 1);
%             this.addParam('Angle offset', 0, 1);
            
            % This class implements the ACoreFun interface!
            this.f = this;
        end
        
        function plot(this, model, t, y)
            plot(y(2,:),y(1,:),'r');
            plot@models.BaseDynSystem(this,model,t,y);
        end
        
        function fx = evaluate(~, x, ~, mu)
            % Implements ACoreFun.evaluate
            a = mu(1);
            b = a+mu(2);
            %b = a+mu(2)+sin(t)/2;
            A = [cos(a) -sin(b); 
                 sin(a) cos(b)];
            fx = A*x;
        end
        
        function varargout = project(varargin)%#ok
            % Implements general.AProjectable
            % no projection for this model as there are only two
            % dimensions.
            error('not implemented');
        end
    end
    
end
