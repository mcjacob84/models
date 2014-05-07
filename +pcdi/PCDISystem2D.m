classdef PCDISystem2D < models.pcdi.BasePCDISystem
    % PCDISystem2D The inhibited version programmed cell death model for 2D
    % geometry.
    %
    % Geometry layout is with zero in the top left corner and x+ being
    % downwards and y+ towards right (according to matlab matrix indexing)
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
        % The kernel expansion used to construct the diffusion coefficient
        % distribution.
        %
        % @type kernels.KernelExpansion
        Kexp;
        
        % Gamma values to use as kernel width for diffusivity
%         Gammas = .08;

        Gammas = [.01 .02 .03 .04 .05];
        
        DiscreteCXMU;
        CurCXMU;
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
            
            nc = 5;
            k.Centers.xi = [.25 1   1.2 .8 .4
                            .45 .25 .65 .9 .8];
            k.Ma = s.rand(1,nc)*1e5;
            this.Kexp = k;
            
            % Spatial area (unscaled!)
            this.Omega = [0 1.5; 0 1] * this.Model.L;
            
            % Scaled!
            this.h = .5 * this.Model.L;
        end
        
        function setConfig(this, mu, inputidx)
            setConfig@models.pcdi.BasePCDISystem(this, mu, inputidx)
            %% Precompute diffusivity c(x,mu) on current grid
            cxmu = this.DiscreteCXMU.compose(0,mu);
            % Use same linear indexing as other quantities on grid
            this.CurCXMU = cxmu(:);
        end
        
        function varargout = plot(this, model, t, y, varargin)
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
        
        function varargout = plotState(this, model, t, y, mode)
            % 
            % Mode 1: Plot line through cell center over time
            % Mode 2: Use an imagesc plot from above with colorbars aside
            % (with pauses)
            % Mode 3: Plot a surface with augmented concentration
            % differences from the mean value (with pauses)
            
            if nargin < 5
%                 mode = [2 3];
                mode = 2;
            end

            % Plot all fields into the same figure
            
            varargout = {};
            if any(mode == 1)
                this.plot1DState(model, t, y, getPM);
            end
            if any(mode == 2) || any(mode == 3)
                pm1 = [];
                if any(mode == 2) 
                    pm1 = getPM;
                end
                pm2 = [];
                if any(mode == 3) 
                    pm2 = getPM;
                end
                this.plot2DState(t, y, pm1, pm2);
            end              
            for k=1:length(varargout)
                varargout{k}.done;
            end
            
            function pm = getPM
                if model.WithInhibitors
                    np = 4;
                else
                    np = 2;
                end
                pm = PlotManager(false,2,np);
                % Use max 2 figures
                pm.MaxFigures = 1;
                pm.LeaveOpen = true;
                varargout{end+1} = pm;
            end
        end
        
        function pm = plotDiffusionCoeff(this, mu, pm, allparts)
            if nargin < 4
                allparts = false;
            end
            if nargin < 3 || isempty(pm)
                if allparts
                    pm = PlotManager(false,2,2);
                else
                    pm = PlotManager;
                end
                if nargout == 0
                    pm.LeaveOpen = true;
                end
                if nargin < 2
                    mu = linspace(0,1,10);
                end
            end
            o = this.Omega / this.Model.L;
            fineness = 1;
            [x,y] = meshgrid(0:this.hs/fineness:o(1,2),0:this.hs/fineness:o(2,2));
            apm = this.getAffParamCX(x,y,'mu');
            if allparts
                for k=1:length(this.Gammas)
                    h = pm.nextPlot(sprintf('diff_coeff_%d',k),...
                        sprintf('Spatial diffusion coefficient base c_%d(x) for \\gamma=%g',...
                        k,this.Gammas(k)),'x','y');
                    surf(h,x,y,log10(apm.getMatrix(k)),'EdgeColor','interp','FaceColor','interp');
                    zlabel('Log-scaled diffusivity');
                end
            end
            for k = 1:length(mu)
                h = pm.nextPlot(sprintf('diff_coeff_mu_%d',k),...
                    sprintf('Spatial diffusion coefficient c(x) for mu(5)=%g',...
                    mu(k)),'x','y');
                surf(h,x,y,log10(apm.compose(0,mu(k))),'EdgeColor','interp','FaceColor','interp');
                zlabel('Log-scaled diffusivity');
            end
            if nargin < 3
                pm.done;
            end
        end
    end
    
    methods(Access=protected)
        function newSysDimension(this)
            m = prod(this.Dims);
            
            %% Compute c(x,mu) on current grid
            o = this.Omega / this.Model.L;
            [x,y] = meshgrid(0:this.hs:o(1,2),0:this.hs:o(2,2));
            this.DiscreteCXMU = this.getAffParamCX(x,y,'mu(5)');
            
            %% Initial conditions
            x0 = dscomponents.AffineInitialValue;
            if this.Model.WithInhibitors
                x0part = zeros(8*m,1);
            else
                x0part = zeros(4*m,1);
            end
            for k = 1:this.DiscreteCXMU.N
                M = this.DiscreteCXMU.getMatrix(k);
                x0part(1:2*m) = [M(:); M(:)]*5e-13;
                x0part(2*m+1:4*m) = [M(:); M(:)]*3e-8;
                if this.Model.WithInhibitors
                    x0part(4*m+1:6*m) = [M(:); M(:)]*3e-8;
                    x0part(6*m+1:end) = [M(:); M(:)]*5e-13;
                end
                x0.addMatrix(this.DiscreteCXMU.funStr{k},x0part);
            end
            this.x0 = x0;
            
            %% Diffusion part
            this.A = this.assembleA;
            
            %% Output extraction
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
            % Pick the Caspase-3 concentration
            ca3 = m+1:2*m;
            C(ca3(sel)) = 1/length(sel);
            this.C = dscomponents.LinearOutputConv(C);
        end
    end
    
    methods(Access=private)
        
        function res = assembleA(this)
            
            D = this.Diff;
            
            res = dscomponents.AffLinCoreFun(this);
            res.TimeDependent = false;
            
            % Add constant diffusion on [0, .2] interval
            A = MatUtils.laplacemat(this.hs, this.Dims(1), this.Dims(2));
            
            % Use coefficient functions that split 1 into n linear
            % functions (equally spaced)
            lso = LinearSplitOfOne('mu(5)',length(this.Gammas));
            add(lso.getFunStr(0),A);
            
            % Then use all specified gamma widths and compute affine-linear
            % interpolation
            for k=1:length(this.Gammas)
                % Set gamma (will influence the functions diffusionCoeff
                % etc)
                this.Kexp.Kernel.Gamma = this.Gammas(k);
                % Compile matrix
                A = this.assembleASpatialDiff;
                add(lso.getFunStr(k),A);
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
        
        function apm = getAffParamCX(this, x, y, muarg)
            % Computes the affine-linear diffusion coefficients on the
            % spatial grid (of the cell) given by x,y
            c = zeros(size(x,1),size(x,2),length(this.Gammas));
            apm = general.AffParamMatrix;
            lso = LinearSplitOfOne(muarg,length(this.Gammas));
            apm.addMatrix(lso.getFunStr(0),ones(size(x)));
            for k=1:length(this.Gammas)
                this.Kexp.Kernel.Gamma = this.Gammas(k);
                hlp = this.diffusionCoeff([x(:) y(:)]');
                c(:,:,k) = reshape(hlp,size(x,1),[]);
                apm.addMatrix(lso.getFunStr(k),c(:,:,k));
            end
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

        function plot2DState(this, t, v, pm1, pm2)
            % Performs a plot for this model's results.
            %
            % Parameters:
            % t: The times `t_0,\ldots,t_N` as row vector @type rowvec
            % v: The system's caspase concentrations (with no output
            % projection!) @type matrix
            
            m = prod(this.Dims);
            pos = reshape(1:8*m,[],8)';

            %% Plot timesteps
            a = this.Omega;
            x = a(1,1):this.h:a(1,2);
            y = a(2,1):this.h:a(2,2);
            [X,Y] = meshgrid(x,y);
            X = X';
            Y = Y';
            
            step = round(length(t)/40);
            for idx=1:step:length(t)
                %% Mode 2 plots
                if ~isempty(pm1)
                    imagescplot(1);
                    imagescplot(2);
                    imagescplot(3);
                    imagescplot(4);
                    if this.Model.WithInhibitors
                        imagescplot(5);
                        imagescplot(6);
                        imagescplot(7);
                        imagescplot(8);
                    end
                    pm1.done;
                end
                set(gcf,'Name',sprintf('Concentration top view at t=%f',t(idx)));
                %% Mode 3 plots
                if ~isempty(pm2)
                    surfplot(1);
                    surfplot(2);
                    surfplot(3);
                    surfplot(4);
                    if this.Model.WithInhibitors
                        surfplot(5);
                        surfplot(6);
                        surfplot(7);
                        surfplot(8);
                    end
                    pm2.done;
                    rotate3d on;
                end
                set(gcf,'Name',sprintf('Augmented concentration differences at t=%f',t(idx)));
                if idx ~= length(t)
                    pause(.1);
                    if any(~ishandle(h))
                        return;
                    end
                end
            end
            
            function imagescplot(dim)
                V = reshape(v(pos(dim,:),idx),this.Dims(1),[])';
                h = pm1.nextPlot(sprintf('PCDI2D_plot2D_%d',this.Tags{dim}),...
                    sprintf('%s concentrations', this.Labels{dim}),'x','y');
                imagesc(x,y,V,'Parent',h);
                set(h,'DataAspectRatio',[1 1 1]);
                axis(h,'xy');
                colorbar('peer',h,'SouthOutside');
            end
            
            function surfplot(dim)
                % Extract the data for the current cell                
                V = reshape(v(pos(dim,:),idx),this.Dims(1),[]);
                
                % Augment the differences
                mv = mean(V(:));
                diff = V - mv;
                diffnorm = norm(diff);
                V = mv + diff/diffnorm;

                h = pm2.nextPlot(sprintf('PCDI2D_plot2D_%d',this.Tags{dim}),...
                    sprintf('%s concentration difference\n Mean %g, augmented by %g',...
                    this.Labels{dim},mv,diffnorm),'x','y');
                %delete(get(h,'Children'));
                surf(h,X,Y,V,'FaceColor','interp','EdgeColor','none');
                ar = get(h,'DataAspectRatio');
                set(h,'DataAspectRatio',[1 1 ar(3)]);
                axis(h,'tight');
            end
        end
    end
    
end

