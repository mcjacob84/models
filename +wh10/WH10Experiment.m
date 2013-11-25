classdef WH10Experiment < models.BaseFullModel
    % Numerical experiments class for Paper WH10
    %
    % Current version works with KerMor 0.4
    
    properties(SetObservable)
        % The system's dimension
        %
        % @propclass{experimental}
        dim;
    end
    
    methods
        
        function this = WH10Experiment(dims)
            this.registerProps('dim');
            
            this.dim = dims;
            
            this.Sampler = [];%sampling.GridSampler;
            
            % This class implements a fake Approx subclass to allow access
            % to the this.Ma property for the error estimator.
            this.Approx = [];
            
            %s = solvers.MLWrapper(@ode45);
            s = solvers.ExplEuler;
            s.MaxStep = [];
            %s = solvers.Heun;
            this.ODESolver = s;
            
            %% System settings
            this.System = models.wh10.WH10System(this);
        end
    end    
    
    methods(Static)
        
        function [d, rmodels] = Experiment1(dim, dir)
            % Experiment with manual reduction and rotation as error source
            %
            % Creates the images also used in the Dissertation (Section
            % 5.1.5, Figure 5.2)
            if nargin < 2
                dir = KerMor.App.TempDirectory;
                if nargin < 1
                    dim = 240000;
                end
            end
            m = models.wh10.WH10Experiment(dim);
            
            m.T = 20;
            m.dt = 0.05;
            
            f = m.System.f;
            
            f.x0 = dscomponents.ConstInitialValue(ones(dim,1));
            f.Ma = repmat(-exp(-f.Centers.xi(1,:)/15),dim,1);
            
            %m.System.Inputs{1} = @(t).4*exp(-t/4); % alter input
            m.System.Inputs{1} = @(t).04*sin(t/3);
            m.System.Inputs{2} = @(t)exp(-abs(10-t/2)); % Gutes Beispiel!
            m.System.Inputs{3} = @(t).5*exp(-(12-t).^2); % Sehr interessante dynamik
            m.System.Inputs{4} = @(t)(t>10)*.1; % okay, aber nicht so spannend
            m.System.B = dscomponents.LinearInputConv(ones(m.dim,1));
            
            %% No training necessary!
            m.TrainingInputs = [];
            
            m.System.C = dscomponents.LinearOutputConv(ones(1,m.dim)/m.dim);
            
            V = ones(dim,1)/sqrt(dim);
            s = spacereduction.ManualReduction(V,V);
            s = spacereduction.RotationDecorator(s);
            s.Dims = 100;
            m.SpaceReducer = s;
            
            d = EstimatorAnalyzer;
            d.EstimatorIterations = [1 2 5];
            d.EstimatorVersions = [1 1 0 0 1 0 0 1 1];
            
            degs = [.0005 .005 .05 .5];
            rmodels = {};
            
            pm = PlotManager;
            pm.AutoTickMarks = false;
            pm.UseFileTypeFolders = false;
            pm.NoTitlesOnSave = true;
            
            %% Create estimation plots
            for didx = 1:length(degs)
                s.Degree = degs(didx);
                d.setModel(m);
                
                rmodels{didx} = d.ReducedModel;%#ok<*AGROW>

                for idx=1:m.System.InputCount
                    m.Name = sprintf('\\theta=%.4f',s.Degree);
                    [errs, relerrs, ctimes] = d.compute([],idx);
                    d.createPlots(errs, relerrs, ctimes, pm);
                    
                    md = d.ModelData(end);
                    % 1 = full, 7 = LSLE TD
                    absmax = 30*abs(md.ErrT(7)-md.ErrT(1));
                    relmax = 30*abs(md.RelErrT(7)-md.RelErrT(1));
                    
                    axis(gca(1),[0 m.T md.MinErr*.9 min(1e4,absmax)]);
                    axis(gca(2),[0 m.T md.MinRelErr*.9 min(1,relmax)]);
                    set(gca(3),'XScale','log');
                    
                    pm.FilePrefix = sprintf('WH10_in%d_deg%f',idx,s.Degree);
                    pm.savePlots(dir,'Format','eps','Close',true);
                end
            end
            
            %% Create output error plot for u3 and theta=0.05
            ma = ModelAnalyzer(rmodels{3});
            ma.SingleFigures = true;
            d.ReducedModel.ErrorEstimator = d.Est(7).Estimator;
            ma.analyzeError([],3,pm);
            pm.FilePrefix = 'WH10_in3_deg05';
            pm.savePlots(dir,'Format','eps','Close',true,'Selection',3);
            
            %% Table stuff
            range = 1:4:13;
            sort = [range (range)+1 (range)+2 (range)+3];
            d.createStatsTables(sort);
        end
        
        function getExp1SysPlots(rmodels)
            ma = ModelAnalyzer;
            ma.SingleFigures = true;
            
            % Pick 0.05 model
            r = rmodels(3);
            h = ma.analyze(r,[],3);
            title(h(2),'');
            Utils.saveFigure(f,'Sys_U3','eps');
        end
        
        function d = Experiment4(dim)
            
            if nargin < 1
                dim = 240000;
            end
            m = WH10Experiment(dim);
            
            m.T = 20;
            m.dt = 0.05;
            
            s = spacereduction.PODReducer;
            s.Mode = 'abs';
            s.Value = 4;
            s.UseSVDS = dim > 10000;
            s = spacereduction.RotationDecorator(s);
            s.Degree = .01;
            s.Dims = 5;
            m.SpaceReducer = s;
            
            f = m.System.f;
            f.x0 = @(mu)ones(dim,1);
            f.Ma = repmat(-exp(-f.Centers.xi(1,:)/15),dim,1);
            
            m.System.Inputs{1} = @(t).2*sin(t);
            %B = rnd.rand(m.dim,1);
            B = linspace(0,1,m.dim)';
            m.System.B = dscomponents.LinearInputConv(B);
            
            m.TrainingInputs = 1;
            
            m.Name = 'Sinus input, mean output mapping, larger error';
            m.System.C = dscomponents.LinearOutputConv(ones(1,m.dim)/m.dim);
            
            d = EstimatorAnalyzer;
            d.EstimatorIterations = [1 2 5];
            d.EstimatorVersions = [1 0 0 1 0 0 1 1];
            d.SingleFigures = true;
            d.setModel(m);
            d.start([],1);
            
            f = d.Figures{1};
            a = gca(f);
            axis(a,[0 m.T 0 10]); 
            set(legend(a),'Location','NorthEast');
            Utils.saveFigure(f,'WH10_e1_errors','fig');
            title(a,'');
            Utils.saveFigure(f,'WH10_e1_errors');
            close(f);
            
            f = d.Figures{2}; a = gca(f);
            axis(a,[0 m.T 5e-5 1]);
            set(legend(a),'Location','NorthEast');
            Utils.saveFigure(f,'WH10_e1_relerrors','fig');
            title(a,'');
            Utils.saveFigure(f,'WH10_e1_relerrors');
            close(f);
            
            f = d.Figures{3}; a = gca(f);
            Utils.saveFigure(f,'WH10_e1_ctimes','fig');
            title(a,'');
            Utils.saveFigure(f,'WH10_e1_ctimes');
            close(f);
        end
        
        function d = Experiment2(dim)
            % Experiment with POD reduction and training with all 4 inputs!
            if nargin < 1
                dim = 240000;
            end
            m = WH10Experiment(dim);
            
            m.T = 20;
            m.dt = 0.05;
            
            f = m.System.f;
            f.x0 = @(mu)ones(dim,1);
            f.Ma = repmat(-exp(-f.Centers.xi(1,:)/15),dim,1);
            
            %m.System.Inputs{1} = @(t).4*exp(-t/4); % alter input
            m.System.Inputs{1} = @(t).04*sin(t/3);
            m.System.Inputs{2} = @(t)exp(-abs(10-t/2)); % Gutes Beispiel!
            m.System.Inputs{3} = @(t)exp(-(10-t).^2); % Sehr interessante dynamik
            m.System.Inputs{4} = @(t)(t>10)*.1; % okay, aber nicht so spannend
            m.System.B = dscomponents.LinearInputConv(ones(m.dim,1));
            
            m.TrainingInputs = 1:4;
            
            m.System.C = dscomponents.LinearOutputConv(ones(1,m.dim)/m.dim);
            
            s = spacereduction.PODReducer;
            s.Value = 4;
            s.Mode = 'abs';
            s.UseSVDS = dim > 10000;
            m.SpaceReducer = s;
            
            d = EstimatorAnalyzer;
            d.EstimatorIterations = [1 2 5];
            d.EstimatorVersions = [1 1 0 0 1 0 0 1 1];
            d.SingleFigures = true;
            d.setModel(m);

            for idx=1:m.System.InputCount
                m.Name = 'POD4in';
                d.start([],idx);

                md = d.ModelData(end);
                % 1 = full, 7 = LSLE TD
                absmax = 30*abs(md.ErrT(7)-md.ErrT(1));
                relmax = 30*abs(md.RelErrT(7)-md.RelErrT(1));

                f = d.Figures{1};
                a = gca(f);
                axis(a,[0 m.T md.ErrT(1)*.9 min(1e4,absmax)]); 
                set(legend(a),'Location','NorthEast');
                Utils.saveFigure(f,sprintf('WH10_POD4_in%d_errors',idx),'fig');
                title(a,'');
                Utils.saveFigure(f,sprintf('WH10_POD4_in%d_errors',idx));

                f = d.Figures{2}; a = gca(f);
                axis(a,[0 m.T md.RelErrT(1)*.9 min(1,relmax)]);
                set(legend(a),'Location','NorthEast');
                Utils.saveFigure(f,sprintf('WH10_POD4_in%d_relerr',idx),'fig');
                title(a,'');
                Utils.saveFigure(f,sprintf('WH10_POD4_in%d_relerr',idx));

                f = d.Figures{3}; a = gca(f);
                Utils.saveFigure(f,sprintf('WH10_POD4_in%d_ctimes',idx),'fig');
                title(a,'');
                Utils.saveFigure(f,sprintf('WH10_POD4_in%d_ctimes',idx));

                close([d.Figures{:}]);
            end
            
            d.createStatsTables;
        end
        
        function d = Experiment3(dim)
            % Experiment with truly random component-wise functions.
            rnd = RandStream('mt19937ar','Seed',2564); %2564
            if nargin < 1
                dim = 240000;
            end
            m = WH10Experiment(dim);
            
            m.T = 20;
            m.dt = 0.05;
            
             s = spacereduction.PODReducer;
            s.Mode = 'abs';
            s.Value = 4;
            s.UseSVDS = dim > 10000;
            s = spacereduction.RotationDecorator(s);
            s.Degree = .01;
            s.Dims = 5;
            m.SpaceReducer = s;
            
            f = m.System.f;
            f.x0 = @(mu)ones(dim,1);
            n = size(f.Centers.xi,2);
            f.Ma = -exp(-rand(dim,n));
            
            m.System.Inputs{1} = @(t).1*sin(t);
            m.System.B = dscomponents.LinearInputConv(rnd.rand(m.dim,1));
            
            this.TrainingInputs = 1;
            
            m.Name = 'Sinus input, mean output mapping, larger error';
            m.System.C = dscomponents.LinearOutputConv(ones(1,m.dim)/m.dim);
            
            d = EstimatorAnalyzer;
            d.EstimatorIterations = [1 2 5];
            d.EstimatorVersions = [1 0 0 1 0 0 1 1];
            d.SingleFigures = true;
            d.setModel(m);
            d.start([],1);
            
%             f = d.Figures{1};
%             a = gca(f);
%             axis(a,[0 m.T 1e-5 1e4]); 
%             set(legend(a),'Location','NorthEast');
%             Utils.saveFigure(f,'WH10_e3_errors','fig');
%             title(a,'');
%             Utils.saveFigure(f,'WH10_e3_errors');
%             close(f);
%             
%             f = d.Figures{2}; a = gca(f);
%             axis(a,[0 m.T 1e-5 1]);
%             set(legend(a),'Location','SouthWest');
%             Utils.saveFigure(f,'WH10_e3_relerrors','fig');
%             title(a,'');
%             Utils.saveFigure(f,'WH10_e3_relerrors');
%             close(f);
%             
%             f = d.Figures{3}; a = gca(f);
%             Utils.saveFigure(f,'WH10_e3_ctimes','fig');
%             title(a,'');
%             Utils.saveFigure(f,'WH10_e3_ctimes');
%             close(f);
        end
         
    end
end

