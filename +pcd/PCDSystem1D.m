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
    
    properties(SetAccess=private, Dependent)
        % Spatial area
        Range;
    end
    
    properties(Access=private)
        % Spatial area
        SRange;
    end
    
    properties(SetAccess=private)
        % System's dimension
        dim;
        
        fRange;
    end
    
    methods
        function this = PCDSystem1D(model)
            this = this@models.pcd.BasePCDSystem(model);

            % Spatial resolution
            elems = 25;
            this.Range = [0 this.Model.L];
            this.h = (this.Range(2)-this.Range(1))/(elems-1);
           
            % Add input param (is getting inserted after the BasePCDSystem
            % condtructor params, so number 9!)
            this.addParam('U', [0.001, 0.1], 12);
            
            %% Initial value
            m = this.dim;
            x0 = zeros(4*m,1);
            x0(m) = this.h;
            
            % Initial 
            %x0(1:15) = 2e-8;
            %x0(m+(1:15)) = 2e-8;
            
            %x0(2*m+(16:m)) = 2e-8;
            %x0(3*m+(16:m)) = 2e-8;
            
            % Initial procasp-concentrations
            %x0(2*m+1:end) = sin(2*pi*(1:2*m)/(2*m))*.01;
            
            %x0(2*m+1:end) = .3;
            %[X,Y] = meshgrid(1:this.dim2,1:this.dim1);
            %s = sin(X * pi/this.dim1) .* exp(-Y/4)*.5;
            %x0(2*m+1:3*m) = s(:);
            this.x0 = dscomponents.ConstInitialValue(x0);
            
            % Extracts the caspase-3 concentrations from the result
            m = this.dim;
            C = zeros(m,4*m);
            C(:,m+1:2*m) = diag(ones(m,1));
            %this.C = dscomponents.LinearInputConv(sparse(C));

            % Set core function
            this.f = models.pcd.CoreFun1D(this);
        end
        
        function set.Range(this, value)
            this.SRange = value/this.Model.L;
            this.fRange = value;
        end
        
        function r = get.Range(this)
            r = this.fRange;
        end
        
        function plot(this, model, t, y)
            % Performs a plot for this model's results.
            %
            figure;
            X = t;
            Y = this.Range(1):this.h:this.Range(2);
            m = model.System.dim;
            doplot(y(1:m,:),'Caspase-8 (x_a)',1);
            doplot(y(m+1:2*m,:),'Caspase-3 (y_a)',2);
            doplot(y(2*m+1:3*m,:),'Pro-Caspase-8 (x_i)',3);
            doplot(y(3*m+1:end,:),'Pro-Caspase-3 (y_i)',4);
            
            function doplot(y,thetitle,pnr)
                subplot(2,2,pnr);
                %mesh(X,Y,y);
                surf(X,Y,y,'EdgeColor','none');
                xlabel('Time [s]');
                ylabel(sprintf('0 to %d: Cell core to hull [m]',this.Range(2)));
                grid off;
                axis tight;
                title(sprintf('Model "%s"\n%s concentrations', model.Name, thetitle));
            end
        end
    end
    
    methods(Access=protected)                
        function updateDims(this)
            % Should not matter which version as rescaling applies to all values
            %m = length(this.SRange(1):this.hs:this.SRange(2));
            m = length(this.Range(1):this.h:this.Range(2));
            this.dim = m;
            % Set state scaling
            ss = zeros(4*m,1);
            ss(1:m) = this.xa0;
            ss(m+1:2*m) = this.ya0;
            ss(2*m+1:3*m) = this.xi0;
            ss(3*m+1:end) = this.yi0;
            this.StateScaling = ss;
        end
    end
    
%     methods(Static,Access=protected)
%         function obj = loadobj(s)
%             obj = models.pcd.PCDSystem1D;
%             obj = loadobj@models.pcd.BasePCDSystem(s, obj);
%             ALoadable.loadProps(mfilename('class'), obj, s);
%         end
%     end
    
end

