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
    
    methods
        function this = PCDSystem3D(model)
            this = this@models.pcd.BasePCDSystem(model);
            
            % Set core function
            this.f = models.pcd.CoreFun3D(this);
            
            % Spatial area [X, Y, Z]
            this.Omega = [0 1; 0 1; 0 1] * this.Model.L;
            this.h = this.Model.L/6;
            
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
        end
      
        function plot(this, model, t, v)
            % Performs a plot for this model's results.
            %
            % Parameters:
            % t: The times `t_0,\ldots,t_N` as row vector @type rowvec
            % v: The system's caspase concentrations (with no output
            % projection!) @type matrix
            
            autocols = true;
            plotback = false;
            
            m = prod(this.Dims);
            xa = v(1:m,:);
            ya = v(m+1:2*m,:);
            xi = v(2*m+1:3*m,:);
            yi = v(3*m+1:end,:);
            b = [min(xa(:)) max(xa(:)); min(ya(:)) max(ya(:));...
                 min(xi(:)) max(xi(:)); min(yi(:)) max(yi(:))];
            %% Prepare figures
            f1 = figure(1); rotate3d on;
            f2 = figure(2); rotate3d on;
            hlpf = figure('Visible','off','MenuBar','none','ToolBar','none');
            hlpax = newplot(hlpf);
            ax = [];
            caps = {'Caspase-8','Caspase-3','Pro-Caspase-8','Pro-Caspase-3'};
            for p = 1:4
                figure(f1);
                a = subplot(2,2,p);
                cla(a);
                set(a,'Box','on','GridLineStyle','none');
                if ~autocols
                    set(a,'CLimMode','manual','CLim',b(p,:));
                end
                view(a,-22,31);
                axis ij;
                axis(a,reshape(this.Omega',1,[]));
                xlabel(a,'Left to right');
                ylabel(a,'Top to bottom');
                zlabel(a,'Front to back');
                title(a,sprintf('Model "%s", %s concentrations', model.Name, caps{p}));
                colorbar('peer',a);
                ax(end+1) = a;%#ok
                
                figure(f2);
                a = subplot(2,2,p);
                set(a,'Box','on','GridLineStyle','none');
                if ~autocols
                    set(a,'CLimMode','manual','CLim',b(p,:));
                end
                view(a,-22,31);
                axis ij;
                axis(a,reshape(this.Omega',1,[]));
                xlabel(a,'Left to right');
                ylabel(a,'Top to bottom');
                zlabel(a,'Front to back');
                title(a,sprintf('Model "%s", %s concentrations', model.Name, caps{p}));
                colorbar('peer',a);
                ax(end+1) = a;%#ok
            end
            
            %% Plot timesteps
            a = this.Omega;
            x = a(1,1):this.h:a(1,2);
            y = a(2,1):this.h:a(2,2);
            z = a(3,1):this.h:a(3,2);
            [X,Y,Z] = meshgrid(x,y,z);
            
            step = round(length(t)/40);
            for idx=1:step:length(t)
                h1 = doplot(xa(:,idx),ax(1:2));
                h2 = doplot(ya(:,idx),ax(3:4));
                h3 = doplot(xi(:,idx),ax(5:6));
                h4 = doplot(yi(:,idx),ax(7:8));
                set(f1,'Name',sprintf('Plot at t=%f',t(idx)));
                set(f2,'Name',sprintf('Plot at t=%f',t(idx)));
                if idx ~= length(t)
                    pause(1);
                    delete([h1; h2; h3; h4]);
                end
            end
            
            close(hlpf);
            
            function h = doplot(zd, ax)        
                V = reshape(zd,this.Dims(1),this.Dims(2),[]);
                %permute(V,[1 3 2]);
                
                hc = contourslice(ax(1),X,Y,Z,V,[],0:.2:1,[]); % 'EdgeColor','none' 
                
                % bugfix: constant nonzero values cause slice to crash when
                % setting the clim property.
                if V(1) ~= 0 && all(V(1) == V(:))
                    V(1) = 1.0001*V(1);
                end
                hs = slice(hlpax,X,Y,Z,V,.5,.5,.5);
                set(hs,'Parent',ax(2),'FaceColor','interp','EdgeColor','none');
                if plotback
                    hs2 = slice(hlpax,X,Y,Z,V,a(1,2),a(2,1),a(3,1));
                    set(hs2,'Parent',ax(2),'FaceColor','interp','EdgeColor','none');
                    hs = [hs; hs2];
                end
                h = [hc; hs];
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
            
%             p = .1; % 10% of each dimensions span, centered in geometry.
%             d = this.dim1;
%             d1idx = find(abs((1:d) - d/2) <= d/2 * p);
%             d = this.dim2;
%             d2idx = find(abs((1:d) - d/2) <= d/2 * p);
%             d = this.dim3;
%             d3idx = find(abs((1:d) - d/2) <= d/2 * p);
%             [d1,d2,d3] = meshgrid(d1idx,d2idx,d3idx);
%             sel = reshape(sub2ind([this.dim1,this.dim2,this.dim3],d1,d2,d3),1,[]);
%             C = zeros(1,4*m);
%             ca3 = m+1:2*m;
%             C(ca3(sel)) = 1/length(sel);
%             this.C = dscomponents.LinearOutputConv(C);
        end
    end
    
    methods(Static)
        function pcd3d_MMACWKA_20logPar_SR(T, dt)
            % Configuration to try and find the good old approximation with
            % little approximation error
            m = models.pcd.PCDModel(3);
            
            if nargin == 1
                dt = T/2000;
            elseif nargin == 0
                T = 8000;
                dt = 1.5;
            end
            m.T = T; %[s]
            m.dt = dt; %[s]
            m.System.h = m.L/6;
            o = solvers.ode.MLode15i;
            o.AbsTol = 1e-6;
            o.RelTol = 1e-5;
            o.MaxStep = [];
            m.ODESolver = o;
            
            s = sampling.GridSampler;
            s.Spacing = 'log';
            m.Sampler = s;
            
            a = m.Approx;
%             s = approx.selection.DefaultSelector;
            %s = approx.selection.LinspaceSelector;
            s = approx.selection.TimeSelector;
            s.Size = 25000;
            a.TrainDataSelector = s;
            aa = approx.algorithms.MinMaxAdaptiveCWKA;
            aa.MaxExpansionSize = 80;
            aa.MaxRelErr = 1e-5;
            aa.MaxAbsErrFactor = 1e-5;
            aa.CheckMaxErrorPercent = .07;
            aa.InitialCenter = 't0';
            a.Algorithm = aa;
            
            m.Data = data.FileModelData(m);
            
            % Zero initial conditions
            %dim = size(m.System.x0.evaluate([]),1);
            %m.System.x0 = dscomponents.ConstInitialValue(zeros(dim,1));
            
%            m.SpaceReducer = [];
            m.SpaceReducer = spacereduction.PODGreedy;
            m.SpaceReducer.Eps = 1e-9;
            
            a = KerMor.App;
            a.Verbose = 2;

            t1 = m.off1_createParamSamples;
            t2 = m.off2_genTrainingData;
            t3 = m.off3_computeReducedSpace;
            t4 = m.off4_genApproximationTrainData;
            
            %factors = general.Utils.createCombinations([5 10 15 20], [1 2 3 4 5],[.1 .01 .001]);
            factors = general.Utils.createCombinations(3, 10, 2);
            n = size(factors,2);
            K = approx.KernelApprox.empty;
            t5 = zeros(1,n);
            for i = 1:n
                fprintf('Starting approximation run %d of %d..\n',i,n);
                aa.NumGammas = factors(1,i);
                aa.MaxGFactor = factors(2,i);
                aa.MinGFactor = factors(3,i);
                t5(i) = m.off5_computeApproximation;
                K(i) = m.Approx.clone;
            end
            times = [t1 t2 t3 t4 t5];%#ok
            
            file = sprintf('pcd3d_MMACWKA_20logPar_SR_T%d_dt%s', T, strrep(sprintf('%f',dt),'.','_'));
            conf = object2str(m);
            save(file, 'm', 'K', 'times', 'factors', 'conf');
        end
    end
end

