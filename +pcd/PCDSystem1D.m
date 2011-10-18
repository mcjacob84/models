classdef PCDSystem1D < models.pcd.BasePCDSystem
    %PCDSYSTEM1D Summary of this class goes here
    %   Detailed explanation goes here
    %
    % @author Daniel Wirtz @date 
    %
    % @change{0,2,dw,2011-03-09} Implemented the model according to the
    % Paper <a
    % href="http://www.simtech.uni-stuttgart.de/publikationen/prints.php?ID=285" 
    % target="_blank">Death wins against life in a spatially extended
    % apoptosis model</a>
    
    methods
        function this = PCDSystem1D(model)
            this = this@models.pcd.BasePCDSystem(model);

            % Set core function
            this.f = models.pcd.CoreFun1D(this);
            
            % Spatial resolution
            elems = 25;
            this.Omega = [0 this.Model.L];
            this.h = (this.Omega(2)-this.Omega(1))/(elems-1);
           
            % Add input param (is getting inserted after the BasePCDSystem
            % condtructor params, so number 9!)
            this.addParam('U', [0.001, 0.1], 12);
        end

        function plot(this, model, t, y)
            % Performs a plot for this model's results.
            %
            figure;
            X = t;
            Y = this.Omega(1):this.h:this.Omega(2);
            m = model.System.Dims(1);
            doplot(y(1:m,:),'Caspase-8 (x_a)',1);
            doplot(y(m+1:2*m,:),'Caspase-3 (y_a)',2);
            doplot(y(2*m+1:3*m,:),'Pro-Caspase-8 (x_i)',3);
            doplot(y(3*m+1:end,:),'Pro-Caspase-3 (y_i)',4);
            
            function doplot(y,thetitle,pnr)
                subplot(2,2,pnr);
                %mesh(X,Y,y);
                surf(X,Y,y,'EdgeColor','none');
                xlabel('Time [s]');
                ylabel(sprintf('%f to %f: Cell core to hull [m]',this.Omega(1),this.Omega(2)));
                grid off;
                axis tight;
                title(sprintf('Model "%s"\n%s concentrations', model.Name, thetitle));
            end
        end
    end
    
    methods(Access=protected)
        function dimsUpdated(this)
            % Assign fitting initial value
            m = prod(this.Dims);
            
            x0 = zeros(4*m,1);
            x0(1:2*m) = 1e-16;
            x0(2*m+1:end) = 1e-9;
            this.x0 = dscomponents.ConstInitialValue(x0);
            
            % Extracts the caspase-3 concentrations from the result
%             C = zeros(m,4*m);
%             C(:,m+1:2*m) = diag(ones(m,1));
            %this.C = dscomponents.LinearInputConv(sparse(C));
        end
    end
end

