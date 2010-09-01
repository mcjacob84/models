classdef PCDSystem1D < models.pcd.BasePCDSystem
    %PCDSYSTEM1D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Spatial area
        Range = [0 1];
    end
    
    properties(SetAccess=private)
        % System's dimension
        dim;
    end
    
    methods
        function this = PCDSystem1D
            this.h = .03;
            
            % Add params
            this.addParam('U', [0, 2], 12);
            
            % Set core function
            this.f = models.pcd.CoreFun1D(this);
            
            % Update simulation constants
            this.updateSimConstants;
        end
        
        function plot(this, model, t, y)
            % Performs a plot for this model's results.
            %
            figure;
            X = t;
            Y = this.Range(1):this.h:this.Range(2);
            mesh(X,Y,y);
            xlabel('Times');
            ylabel(sprintf('0 to %d: Cell core to hull',this.Range(2)));
            grid off;
            axis tight;
            title(sprintf('Plot for output of model "%s"', model.Name));
        end
    end
    
    methods(Access=protected)
        function x0 = initialX(this, mu)%#ok
            m = this.dim;
            x0 = zeros(4*m,1);
            %x0(2*m+1:end) = .3;
            %[X,Y] = meshgrid(1:this.dim2,1:this.dim1);
            %s = sin(X * pi/this.dim1) .* exp(-Y/4)*.5;
            %x0(2*m+1:3*m) = s(:);
        end
        
        function C = getC(this, t, mu)%#ok
            % Extracts the caspase-3 concentrations from the result
            m = this.dim;
            C = zeros(m,4*m);
            C(:,m+1:2*m) = diag(ones(m,1));
            C = sparse(C);
        end
        
        function updateDimSimConstants(this)
            this.dim = length(this.Range(1):this.h:this.Range(2));
        end
    end
    
end

