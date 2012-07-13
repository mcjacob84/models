classdef PCDSystem2D < models.pcd.BasePCDSystem
    %PCDSystem2D The programmed cell death model for 2D geometry.
    % 
    % The first row of Omega denotes the width of the geometry and
    % the second the height.
    %
    % @change{0,3,sa,2011-05-11} Implemented property setter   
    
    methods
        function this = PCDSystem2D(model)
            this = this@models.pcd.BasePCDSystem(model);
            
            % Set core function
            this.f = models.pcd.CoreFun2D(this);
            
            % Spatial area (unscaled!)
            this.Omega = [0 1.5; 0 1] * this.Model.L;
            % Scaled!
            this.h = .5 * this.Model.L;
            %model.dt = model.dt*5;
            
            % Add params
            % Simpler version: homogeneous area & rate
            this.addParam('area', [0, 1], 10);
            this.addParam('rate', [1e-5, 1e-2], 15);
%             rate_min = 1e-4;
%             rate_max = 1e-1;
%             % Param indices 9-12
%             this.addParam('area_top', [0, .9], 3);
%             this.addParam('area_bottom', [0, 0], 1);
%             this.addParam('area_left', [0, .9], 3);
%             this.addParam('area_right', [0, 0], 1);
%             % Param indices 13-16
%             this.addParam('rate_top', [rate_min, rate_max], 3);
%             this.addParam('rate_bottom', [.005, .005], 1);
%             this.addParam('rate_left', [rate_min, rate_max], 3);
%             this.addParam('rate_right', [rate_min, rate_max], 1); 
        end        

        function plot(this, model, t, v)
            % Performs a plot for this model's results.
            %
            % Parameters:
            % t: The times `t_0,\ldots,t_N` as row vector @type rowvec
            % v: The system's caspase concentrations (with no output
            % projection!) @type matrix
            
            autocols = true;
            
            m = prod(this.Dims);
            xa = v(1:m,:);
            ya = v(m+1:2*m,:);
            xi = v(2*m+1:3*m,:);
            yi = v(3*m+1:end,:);
            b = [min(xa(:)) max(xa(:)); min(ya(:)) max(ya(:));...
                 min(xi(:)) max(xi(:)); min(yi(:)) max(yi(:))];
            %% Prepare figures
            f1 = figure(1); rotate3d on;
            hlpf = figure('Visible','off','MenuBar','none','ToolBar','none');
            hlpax = newplot(hlpf);
            axis tight;
            ax = [];
            caps = {'Caspase-8','Caspase-3','Pro-Caspase-8','Pro-Caspase-3'};
            for p = 1:4
                figure(f1);
                a = subplot(2,2,p);
                cla(a);
                %axis tight;
                axis(a,[reshape(this.Omega',1,[]) b(p,:)]);
                ar = get(gca,'DataAspectRatio');
                set(gca,'DataAspectRatio',[ar(1)/2 ar(2:3)]);
                %axis(a fill;
                set(a,'Box','on','GridLineStyle','none');
                if ~autocols
                    set(a,'CLimMode','manual','CLim',b(p,:));
                end
                view(a,-22,31);
                xlabel(a,'Left to right');
                ylabel(a,'Bottom to Top');
                zlabel(a,sprintf('%s concentration',caps{p}));
                title(a,sprintf('Model "%s", %s concentrations', model.Name, caps{p}));
                colorbar('peer',a);
                ax(end+1) = a;%#ok
            end
            
            %% Plot timesteps
            a = this.Omega;
            x = a(1,1):this.h:a(1,2);
            y = a(2,1):this.h:a(2,2);
            [X,Y] = meshgrid(x,y);
            
            step = round(length(t)/40);
            for idx=1:step:length(t)
                h1 = doplot(xa(:,idx),ax(1));
                h2 = doplot(ya(:,idx),ax(2));
                h3 = doplot(xi(:,idx),ax(3));
                h4 = doplot(yi(:,idx),ax(4));
                set(f1,'Name',sprintf('Plot at t=%f',t(idx)));
                if idx ~= length(t)
                    pause;
                    delete([h1; h2; h3; h4]);
                end
            end
            
            close(hlpf);
            
            function hs = doplot(zd, ax)       
                V = reshape(zd,this.Dims(1),[])';
                                
                % bugfix: constant nonzero values cause slice to crash when
                % setting the clim property.
%                 if V(1) ~= 0 && all(V(1) == V(:))
%                     V(1) = 1.0001*V(1);
%                 end
                %hs = mesh(hlpax,X,Y,V);
                hs = surf(hlpax,X,Y,V);
                set(hs,'Parent',ax);
                set(hs,'FaceColor','interp','EdgeColor','none');
            end
        end
    end
    
    methods(Access=protected)            
        function newSysDimension(this)
            m = prod(this.Dims);
            x0 = zeros(4*m,1);
            x0(1:2*m) = 1e-16;
            x0(2*m+1:end) = 1e-9;
            this.x0 = dscomponents.ConstInitialValue(x0);
            
            % Diffusion part
            A = general.MatUtils.laplacemat(this.hs, this.Dims(1), this.Dims(2));
            A = blkdiag(A,this.Diff(1)*A,this.Diff(2)*A,this.Diff(3)*A);
            this.A = dscomponents.LinearCoreFun(A);
            
            % Output extraction
%             p = .1; % 10% of each dimensions span, centered in geometry.
%             d = this.dim1;
%             d1idx = find(abs((1:d) - d/2) <= d/2 * p);
%             d = this.dim2;
%             d2idx = find(abs((1:d) - d/2) <= d/2 * p);
%             [d1,d2] = meshgrid(d1idx,d2idx);
%             sel = reshape(sub2ind([this.dim1,this.dim2],d1,d2),1,[]);
%             C = zeros(1,4*m);
%             ca3 = m+1:2*m;
%             C(ca3(sel)) = 1/length(sel);
%             this.C = dscomponents.LinearOutputConv(C);
            this.C = dscomponents.LinearOutputConv(1);
        end
    end
    
end

