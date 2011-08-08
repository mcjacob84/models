classdef PCDSystem2D < models.pcd.BasePCDSystem
    %PCDSYSTEM2D Summary of this class goes here
    %   Detailed explanation goes here
    %
    % @change{0,3,sa,2011-05-11} Implemented property setter
    
    properties(Constant)
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
            
            this.h = .1;
            
            % Add params
            this.addParam('reac_rate_upper', [0, 0.002], 5);
            this.addParam('reac_rate__right', 0, 1);
            this.addParam('reac_rate_lower', [0, 0.002], 3);
            this.addParam('reac_rate_left', 0, 1);
            
            m = this.dim1*this.dim2;
            x0 = zeros(4*m,1);
            %x0(2*m+1:end) = .3;
            %[X,Y] = meshgrid(1:this.dim2,1:this.dim1);
            %s = sin(X * pi/this.dim1) .* exp(-Y/4)*.5;
            %x0(2*m+1:3*m) = s(:);
            this.x0 = dscomponents.ConstInitialValue(x0);
            
            p = .1; % 10% of each dimensions span, centered in geometry.
            d = this.dim1;
            d1idx = find(abs((1:d) - d/2) <= d/2 * p);
            d = this.dim2;
            d2idx = find(abs((1:d) - d/2) <= d/2 * p);
            [d1,d2] = meshgrid(d1idx,d2idx);
            sel = reshape(sub2ind([this.dim1,this.dim2],d1,d2),1,[]);
            C = zeros(1,4*m);
            ca3 = m+1:2*m;
            C(ca3(sel)) = 1/length(sel);
            this.C = dscomponents.LinearOutputConv(C);
            
            % Set core function
            this.f = models.pcd.CoreFun2D(this);
        end
        
%         function set.Omega(this,value)
%             if ~isa(value, 'double')
%                 error('Value must be a valid double matrix');
%             end
%             this.Omega = value;
%         end
        
        function plot(this, model, t, z)
            % Performs a plot for this model's results.
            %
            fh = figure;
            %[X,Y] = meshgrid(1:d1,1:d2);
            x = this.Omega(1,1):this.h:this.Omega(1,2);
            y = this.Omega(2,1):this.h:this.Omega(2,2);
            m = this.dim1*this.dim2;
            step = round(length(t)/20);
            for idx=1:step:length(t)
                doplot(z(1:m,idx),'Caspase-8',1);
                doplot(z(m+1:2*m,idx),'Caspase-3',2);
                doplot(z(2*m+1:3*m,idx),'Pro-Caspase-8',3);
                doplot(z(3*m+1:end,idx),'Pro-Caspase-3',4);
                set(fh,'name',sprintf('Plot at t=%f',t(idx)));
                pause;
            end
            
            function doplot(zd,thetitle,pnr)
                subplot(2,2,pnr);
%                 ma = max(zd(:));
%                 if ma ~= 0
%                     zd = zd/ma;
%                 end 
                Z = reshape(zd,this.dim1,[]);
                if pnr == 4
                    Z
                    figure(5);
                    imagesc(Z);
                    figure(fh);
                end
                [X,Y] = meshgrid(x,y);
                mesh(X,Y,Z); % 'EdgeColor','none' 
                xlabel(sprintf('%f to %f',this.Omega(1,1),this.Omega(1,2)));
                ylabel(sprintf('%f to %f',this.Omega(2,1),this.Omega(2,2)));
                %grid off;
                axis tight;
                title(sprintf('Model "%s", %s concentrations', model.Name, thetitle));
            end
        end
    end
    
    methods(Access=protected)
        
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
            m = this.dim1*this.dim2;
            ss = zeros(4*m,1);
            ss(1:m) = this.xa0;
            ss(m+1:2*m) = this.ya0;
            ss(2*m+1:3*m) = this.xi0;
            ss(3*m+1:end) = this.yi0;
            this.StateScaling = ss;
        end
    end
    
end

