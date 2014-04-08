classdef PCDISystem2D < models.pcdi.BasePCDISystem
    %PCDISystem2D The inhibited version programmed cell death model for 2D
    %geometry.
    % 
    % @author Daniel Wirtz @date 2013-10-23
    %
    % @new{0,8,dw,2013-10-23} Added this class
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.morepas.org/software/index.html
    % - \c Documentation http://www.morepas.org/software/kermor/index.html
    % - \c License @ref licensing
    
    properties
        % Flag that indicates if the plot command yields a 2D plot or a 1D slice over time
        % plot.
        %
        % @type logical @default false
        Plot2D = false;
        
        % The kernel expansion used to construct the diffusion coefficient
        % distribution.
        %
        % @type kernels.KernelExpansion
        Kexp;
        
        Gammas = [.03 .06 .08 .1];
    end
    
    methods
        function this = PCDISystem2D(model)
            this = this@models.pcdi.BasePCDISystem(model);
            
            % Set core function
            if model.WithInhibitors
                this.f = models.pcdi.InhibitCoreFun2D(this);
            else
                this.f = models.pcdi.CoreFun2D(this);
            end
            
            % Assemble kernel expansion
            s = RandStream('mt19937ar','Seed',2);
            k = kernels.KernelExpansion;
%             ke = kernels.Wendland;
%             ke.d = 2;
%             ke.k = 3;
%             k.Kernel = ke;
            k.Kernel.Gamma = .08;
            
            nc = 5;
            k.Centers.xi = [.25 1   1.2 .8 .4
                            .45 .25 .65 .9 .8];
            k.Ma = s.rand(1,nc)*10000;
            this.Kexp = k;
            
            % Spatial area (unscaled!)
            this.Omega = [0 1.5; 0 1] * this.Model.L;
            
            % Scaled!
            this.h = .5 * this.Model.L;
        end
        
        function varargout = plot(~, model, t, y, varargin)
            if ~isempty(varargin) && isa(varargin{1},'PlotManager')
                pm = varargin{1};
            else
                pm = PlotManager;
                if nargout == 0
                    pm.LeaveOpen = true;
                else
                    varargout{1} = pm;
                end
            end
            h = pm.nextPlot([model.SaveTag '_outputplot'],...
                sprintf('Output plot for model %s',model.Name),'Time','Caspase-3 Concentration');
            %plot(h,t,y,'r','LineWidth',2);
            semilogy(h,t,y,'r','LineWidth',2);
            if isempty(varargin)
                pm.done;
            end
        end
        
        function varargout = plotState(this, model, t, y, varargin)
            if ~isempty(varargin) && isa(varargin{1},'PlotManager')
                pm = varargin{1};
            else
                if model.WithInhibitors
                    np = 4;
                else
                    np = 2;
                end
                pm = PlotManager(false,2,np);
                if nargout == 0
                    pm.LeaveOpen = true;
                else
                    varargout{1} = pm;
                end
            end
            if this.Plot2D
                this.plot2DState(model, t, y, pm);
            else
                this.plot1DState(model, t, y, pm);
            end
            if isempty(varargin)
                pm.done;
            end
        end
        
        function pm = plotDiffusionCoeff(this, pm)
            if nargin < 2
                pm = PlotManager(false,2,2);
                if nargout == 0
                    pm.LeaveOpen = true;
                end
            end
            o = this.Omega / this.Model.L;
            [x,y] = meshgrid(0:this.hs/15:o(1,2),0:this.hs/15:o(2,2));
            for k=1:length(this.Gammas)
                this.Kexp.Kernel.Gamma = this.Gammas(k);
                c = this.diffusionCoeff([x(:) y(:)]');
                c = reshape(c,size(x,1),[]);
                h = pm.nextPlot('diff_coeff','Spatial diffusion coefficient c(x)','x','y');
                surf(h,x,y,c,'EdgeColor','interp');
            end
        end
    end
    
    methods(Access=protected)
        function newSysDimension(this)
            m = prod(this.Dims);
            if this.Model.WithInhibitors
                x0 = zeros(8*m,1);
            else
                x0 = zeros(4*m,1);
            end
            x0(1:2*m) = 5e-13;
            x0(2*m+1:4*m) = 3e-8;
            if this.Model.WithInhibitors
                x0(4*m+1:6*m) = 3e-8;
                x0(6*m+1:end) = 5e-13;
            end
            this.x0 = dscomponents.ConstInitialValue(x0);
            
            % Diffusion part
            this.A = this.assembleA;
            
            % Output extraction
            p = .1; % 10% of each dimensions span, centered in geometry.
            d = this.Dims(1);
            d1idx = find(abs((1:d) - d/2) <= d/2 * p);
            d = this.Dims(2);
            d2idx = find(abs((1:d) - d/2) <= d/2 * p);
            [d1,d2] = meshgrid(d1idx,d2idx);
            sel = reshape(sub2ind(this.Dims,d1,d2),1,[]);
            if this.Model.WithInhibitors
                C = sparse(1,8*m);
            else
                C = sparse(1,4*m);
            end
            ca3 = m+1:2*m;
            C(ca3(sel)) = 1/length(sel);
            this.C = dscomponents.LinearOutputConv(C);
        end
    end
    
    methods(Access=private)
        
        function res = assembleA(this)
            
            D = this.Diff;
            
            res = dscomponents.AffLinCoreFun;
            res.TimeDependent = false;
            
            ilen = 1/length(this.Gammas);
            % Add constant diffusion on [0, .2] interval
            A = MatUtils.laplacemat(this.hs, this.Dims(1), this.Dims(2));
            add(sprintf('(1-mu(5)/%g)*(mu(5)<%g)',ilen,ilen),A);
            
            % Then use all specified gamma widths and compute affine-linear
            % interpolation
            for k=1:length(this.Gammas)
                % Set gamma (will influence the functions diffusionCoeff
                % etc)
                this.Kexp.Kernel.Gamma = this.Gammas(k);
                % Compile matrix
                A = this.assembleASpatialDiff;
                % zero to one
                zto = sprintf('(mu(5)>=%g)*(mu(5)<%g)*((mu(5)-%g)/%g)',(k-1)*ilen,k*ilen,(k-1)*ilen,ilen);
                % one to zero
                otz = sprintf('(mu(5)>=%g)*(mu(5)<%g)*(1-(mu(5)-%g)/%g)',k*ilen,(k+1)*ilen,k*ilen,ilen);
                add([zto '+' otz],A);
            end
            
            function add(str, A)
                if this.Model.WithInhibitors
                    augA = blkdiag(A,D(1)*A,D(2)*A,D(3)*A,...
                        D(4)*A,D(5)*A,D(6)*A,D(7)*A);
                else
                    augA = blkdiag(A,D(1)*A,D(2)*A,D(3)*A);
                end
                res.addMatrix(str, augA);
            end
            
        end
        
        function A = assembleASpatialDiff(this)
            % Get scaled area
            o = this.Omega / this.Model.L;
            [x,y] = meshgrid(0:this.hs:o(1,2),0:this.hs:o(2,2));
            x = x'; y = y';
            if size(x,1) ~= this.Dims(1) || size(x,2) ~= this.Dims(2)
                error('Boo');
            end
            
            %% c(x) * \laplace u part
            % laplace matrix
            A1 = MatUtils.laplacemat(this.hs, this.Dims(1), this.Dims(2));
            % times c(x)
            A1 = bsxfun(@times,A1,this.diffusionCoeff([x(:) y(:)]')');
            
            %% grad c(x) . grad u part
            A2 = MatUtils.divcdivumat(x,y,@this.nablaC);
            
            A = A1 + A2;
        end
        
        function c = diffusionCoeff(this, x)
            c = 1./(1+this.Kexp.evaluate(x));
        end
        
        function nc = nablaC(this, x)
            nc = -1/(1+this.Kexp.evaluate(x))^2 * this.Kexp.getStateJacobian(x);
        end
        
        function plot1DState(this, model, t, y, pm)
            m = prod(this.Dims);
            
            % Select cell center values
            idxmat = zeros(this.Dims);
            idxmat(:) = 1:m;
            sel = idxmat(:,round(this.Dims(2)/2));
            sel = [sel; sel+m; sel+2*m; sel+3*m;...
                sel+4*m; sel+5*m; sel+6*m; sel+7*m];
            m = length(sel)/8;
            y = y(sel,:);
            
            if length(t) > 150
                idx = round(linspace(1,length(t),150));
                t = t(idx);
                y = y(:,idx);
            end
            states = {'dead', 'alive', 'unstable'};
            ss = this.Model.getSteadyStates(this.mu(4))';
            
            X = t;
            Y = (this.Omega(1,1):this.h:this.Omega(1,2))/model.L;
            pos = reshape(1:8*m,[],8)';
            doplot('c8','Caspase-8 (x_a)',1);
            doplot('c3','Caspase-3 (y_a)',2);
            doplot('pc8','Pro-Caspase-8 (x_i)',3);
            doplot('pc3','Pro-Caspase-3 (y_i)',4);
            if this.Model.WithInhibitors
                doplot('iap','IAP (iap)',5);
                doplot('bar','BAR (bar)',6);
                doplot('yb','Caspase-3+IAP (yb)',7);
                doplot('xb','Caspase-8+BAR (xb)',8);
            end
            
            function doplot(tag, thetitle, pnr)
                yl = y(pos(pnr,:),:);
                perc = (yl(end) - ss(2,pnr)) / (ss(1,pnr) - ss(2,pnr));
                if perc < .5
                    id = 2;
                else
                    id = 1;
                end
                    
%                 if any(reldi > .1) || any(reldi < 10)
%                     [~, id] = min(perc);
                    tit = sprintf('%s concentrations\nCell state at T=%.4g: %s\n%.5g (%5.2f%%)', thetitle,...
                    max(t),states{id},yl(end),perc*100);
%                 else
%                     tit = sprintf('%s concentrations\n%s', thetitle,reldistr);
%                 end
                if pm.Single
                    tit = sprintf('Model "%s"\n%s',model.Name,tit);
                end
                h = pm.nextPlot(tag,tit,'Time [s]','Cell slice');
                surf(h,X,Y,yl,'EdgeColor','none');
                zlabel(h,thetitle);
            end
        end

        function plot2DState(this, model, t, v, pm)
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
            rotate3d on;
            hlpf = figure('Visible','off','MenuBar','none','ToolBar','none');
            hlpax = newplot(hlpf);
            axis tight;
            ax = [];
            caps = {'Caspase-8','Caspase-3','Pro-Caspase-8','Pro-Caspase-3'};
            for p = 1:4
                h = pm.nextPlot(sprintf('PCD2D_plot1D_%d',p),...
                    sprintf('Model "%s", %s concentrations', model.Name, caps{p}),...
                    'Left to right','Bottom to Top');
                axis(h,[reshape(this.Omega',1,[]) b(p,:)]);
                ar = get(gca,'DataAspectRatio');
                set(h,'DataAspectRatio',[ar(1)/2 ar(2:3)]);
                %axis(a fill;
                set(h,'Box','on','GridLineStyle','none');
                if ~autocols
                    set(h,'CLimMode','manual','CLim',b(p,:));
                end
                view(h,-22,31);
                zlabel(h,sprintf('%s concentration',caps{p}));
                colorbar('peer',h);
                ax(end+1) = h;%#ok
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
                set(gcf,'Name',sprintf('Plot at t=%f',t(idx)));
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
    
end

