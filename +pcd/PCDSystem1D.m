classdef PCDSystem1D < models.pcd.BasePCDSystem
    %PCDSYSTEM1D Summary of this class goes here
    %   Detailed explanation goes here
    %
    % @author Daniel Wirtz @date 
    %
    % @change{0,3,dw,2011-03-09} Implemented the model according to the
    % Paper <a
    % href="http://www.simtech.uni-stuttgart.de/publikationen/prints.php?ID=285" 
    % target="_blank">Death wins against life in a spatially extended
    % apoptosis model</a>
    
    properties
        % Spatial area
        Range = [0 1];
    end
    
    properties(SetAccess=private)
        % System's dimension
        dim;
    end
    
    methods
        function this = PCDSystem1D(model)
            this = this@models.pcd.BasePCDSystem(model);
            
            this.h = .03;
            
            % Add input param (is getting inserted after the BasePCDSystem
            % condtructor params, so number 9!)
            this.addParam('U', [0, 0.001], 12);
            
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
            m = model.System.dim;
            doplot(y(1:m,:),'Caspase-8',1);
            doplot(y(m+1:2*m,:),'Caspase-3',2);
            doplot(y(2*m+1:3*m,:),'Pro-Caspase-8',3);
            doplot(y(3*m+1:end,:),'Pro-Caspase-3',4);
            
            function doplot(y,thetitle,pnr)
                subplot(2,2,pnr);
                %mesh(X,Y,y);
                surf(X,Y,y,'EdgeColor','none');
                xlabel('Time [s]');
                ylabel(sprintf('0 to %d: Cell core to hull [m]',this.Range(2)));
                grid off;
                axis tight;
                title(sprintf('Model "%s", %s concentrations', model.Name, thetitle));
            end
        end
    end
    
    methods(Access=protected)
        function x0 = initialX(this, mu)%#ok
            m = this.dim;
            x0 = zeros(4*m,1);
            
            % Initial 
            %x0(round(2*m/3):m) = .001;
            
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
        
        function updateDims(this)
            this.dim = length(this.Range(1):this.h:this.Range(2));
        end
    end
    
end

