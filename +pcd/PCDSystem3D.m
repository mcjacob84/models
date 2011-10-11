classdef PCDSystem3D < models.pcd.BasePCDSystem
%PCDSYSTEM3D 3D implementation of the cell apoptosis model
%
% @author Daniel Wirtz @date 2011-10-05
%
% @new{0,5,dw,2011-10-05} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(Constant)
        % Spatial area
        Omega = [0 1; 0 1; 0 1];
    end
    
    properties(SetAccess=private)
        % System's dimension 1 along x
        dim1;
        
        % System's dimension 2 along y
        dim2;
        
        % System's dimension 3 along z
        dim3;
    end
    
    methods
        function this = PCDSystem3D(model)
            this = this@models.pcd.BasePCDSystem(model);
            
            this.h = .1;
            
            % Add params
            rate_max = 1e-4;
            % Param indices 9-14
            this.addParam('area_front', [.3, .8], 3);
            this.addParam('area_back', [0, 0], 2);
            this.addParam('area_left', [.3, .8], 3);
            this.addParam('area_right', [0, 0], 2);
            this.addParam('area_top', [.3, .8], 2);
            this.addParam('area_bottom', [0, 0], 2);
            % Param indices 15-20
            this.addParam('rate_front', [0, rate_max], 4);
            this.addParam('rate_back', [0, rate_max], 1);
            this.addParam('rate_left', [0, rate_max], 4);
            this.addParam('rate_right', [0, rate_max], 1);
            this.addParam('rate_top', [0, rate_max], 4);
            this.addParam('rate_bottom', [0, rate_max], 1);
            
            m = this.dim1*this.dim2*this.dim3;
            x0 = zeros(4*m,1);
            this.x0 = dscomponents.ConstInitialValue(x0);
            
            p = .1; % 10% of each dimensions span, centered in geometry.
            d = this.dim1;
            d1idx = find(abs((1:d) - d/2) <= d/2 * p);
            d = this.dim2;
            d2idx = find(abs((1:d) - d/2) <= d/2 * p);
            d = this.dim3;
            d3idx = find(abs((1:d) - d/2) <= d/2 * p);
            [d1,d2,d3] = meshgrid(d1idx,d2idx,d3idx);
            sel = reshape(sub2ind([this.dim1,this.dim2,this.dim3],d1,d2,d3),1,[]);
            C = zeros(1,4*m);
            ca3 = m+1:2*m;
            C(ca3(sel)) = 1/length(sel);
            this.C = dscomponents.LinearOutputConv(C);
            
            % Set core function
            this.f = models.pcd.CoreFun3D(this);
        end
      
        function plot(this, model, t, v)
            % Performs a plot for this model's results.
            %
            fh = figure;
            %[X,Y] = meshgrid(1:d1,1:d2);
            x = this.Omega(1,1):this.h:this.Omega(1,2);
            y = this.Omega(2,1):this.h:this.Omega(2,2);
            z = this.Omega(3,1):this.h:this.Omega(3,2);
            m = this.dim1*this.dim2*this.dim3;
            step = round(length(t)/20);
            for idx=1:step:length(t)
                doplot(v(1:m,idx),'Caspase-8',1);
                doplot(v(m+1:2*m,idx),'Caspase-3',2);
                doplot(v(2*m+1:3*m,idx),'Pro-Caspase-8',3);
                doplot(v(3*m+1:end,idx),'Pro-Caspase-3',4);
                set(fh,'name',sprintf('Plot at t=%f',t(idx)));
                pause;
            end
            
            function doplot(zd,thetitle,pnr)
                subplot(2,2,pnr); cla;
%                 ma = max(zd(:));
%                 if ma ~= 0
%                     zd = zd/ma;
%                 end 
                V = reshape(zd,this.dim1,this.dim2,[]);
                [X,Y,Z] = meshgrid(x,y,z);
                contourslice(X,Y,Z,V,[],0:.2:1,[]); % 'EdgeColor','none' 
                xlabel(sprintf('%f to %f',this.Omega(1,1),this.Omega(1,2)));
                ylabel(sprintf('%f to %f',this.Omega(2,1),this.Omega(2,2)));
                zlabel(sprintf('%f to %f',this.Omega(3,1),this.Omega(3,2)));
                %grid off;
                axis tight;
                title(sprintf('Model "%s", %s concentrations', model.Name, thetitle));
            end
        end
    end
    
    methods(Access=protected)        
        function updateDims(this)
            this.dim1 = length(this.Omega(1,1):this.h:this.Omega(1,2));
            this.dim2 = length(this.Omega(2,1):this.h:this.Omega(2,2));
            this.dim3 = length(this.Omega(3,1):this.h:this.Omega(3,2));
            m = this.dim1*this.dim2*this.dim3;
            ss = zeros(4*m,1);
            ss(1:m) = this.xa0;
            ss(m+1:2*m) = this.ya0;
            ss(2*m+1:3*m) = this.xi0;
            ss(3*m+1:end) = this.yi0;
            this.StateScaling = ss;
        end
    end
    
end

