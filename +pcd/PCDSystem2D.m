classdef PCDSystem2D < models.pcd.BasePCDSystem
    %PCDSYSTEM2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Spatial area
        Omega = [0 1; 0 1];
    end
    
    properties(SetAccess=private)
        % System's dimension 1
        dim1;
        % System's dimension 2
        dim2;
    end
    
    methods
        function this = PCDSystem2D(model)
            this = this@models.pcd.BasePCDSystem(model);
            
            this.h = .2;
            
            % Add params
            this.addParam('mu1', [0, 2], 5);
            this.addParam('mu2', [0, 2], 3);
            this.addParam('mu3', [0, 1], 3);
            this.addParam('mu4', [0, 1], 3);
            
            % Set core function
            this.f = models.pcd.CoreFun2D(this);
            
            % Update simulation constants
            this.updateSimConstants;
        end
        
        function plot(this, model, t, y)
            % Performs a plot for this model's results.
            %
            figure;
            %[X,Y] = meshgrid(1:d1,1:d2);
            X = this.Omega(1,1):this.h:this.Omega(1,2);
            Y = this.Omega(2,1):this.h:this.Omega(2,2);
            for idx=1:length(t)
                Z = reshape(y(:,idx),this.dim1,[]);
                surf(X,Y,Z);
                grid off;
                %imagesc(Z);
                axis tight;
                title(sprintf('Plot for output of model "%s" at time %1.2f', model.Name, t(idx)));
                pause;
            end
        end
    end
    
    methods(Access=protected)
        function x0 = initialX(this, mu)%#ok
            m = this.dim1*this.dim2;
            x0 = zeros(4*m,1);
            %x0(2*m+1:end) = .3;
            %[X,Y] = meshgrid(1:this.dim2,1:this.dim1);
            %s = sin(X * pi/this.dim1) .* exp(-Y/4)*.5;
            %x0(2*m+1:3*m) = s(:);
        end
        
        function C = getC(this, t, mu)%#ok
            % Extracts the caspase-3 concentrations from the result
            m = this.dim1*this.dim2;
            C = zeros(m,4*m);
            C(:,m+1:2*m) = diag(ones(m,1));
            C = sparse(C);
        end
        
        function updateDims(this)
            this.dim1 = length(this.Omega(1,1):this.h:this.Omega(1,2));
            this.dim2 = length(this.Omega(2,1):this.h:this.Omega(2,2));
        end
    end
    
end

