classdef MathMODExperiment < models.BaseFullModel
    % Numerical experiments class for Paper ICIAM
    %
    % Current version works with KerMor 0.4
    
    properties(SetObservable)
        % The system's dimension
        %
        % @propclass{experimental}
        dim;
    end
    
    methods
        
        function this = MathMODExperiment(dims)
            this.registerProps('dim');
            
            if nargin == 0
                dims = 100;
            end
            this.dim = dims;
            
            this.Name = 'MathMOD 2012 Experiments';
            
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
            this.System = models.mathmod2012.MathMODSystem(this);
            
            this.T = 20;
            this.dt = 0.05;
            this.System.MaxTimestep = this.dt;
            
            this.Sampler = sampling.GridSampler;
            %m.System.addParam('mu3',[3 3],1);
            
            e1 = [ones(dims,1) zeros(dims,1)];
            e2 = [zeros(dims,1) ones(dims,1)];
            
            x0 = dscomponents.AffineInitialValue;
            x0.addMatrix('mu(3)', ones(dims,1));
            %x0.addMatrix('sum(mu)', e1);
            this.System.x0 = x0;
            
            f = this.System.f.Expansion;
            f.Ma = repmat(-exp(-f.Centers.xi(1,:)/15),dims,1);
            
            this.System.Inputs{1} = @(t)[.4*sin(t/3); exp(-(12-t).^2)];
            this.System.Inputs{2} = @(t)[.5*sin(t/2); 4*exp(-7*(12-t).^2)-.5*exp(-(5-t).^2)];
            B = dscomponents.AffLinInputConv;
            B.addMatrix('mu(1)', e1);
            B.addMatrix('1-mu(1)', e2);
            %             B.addMatrix('-mu(2)', e2);
            this.System.B = B;
            
            this.TrainingInputs = 1;
            
            this.System.C = dscomponents.LinearOutputConv(ones(1,dims)/dims);
            
            V = ones(dims,1)/sqrt(dims);
            s = spacereduction.ManualReduction(V,V);
            %             s = spacereduction.PODReducer;
            %             s.Value = 3;
            %             s.Mode = 'abs';
            %             s.UseSVDS = dims > 10000;
            %             this.SpaceReducer = s;
            
            s = spacereduction.RotationDecorator(s);
            s.Dims = 100;
            s.Degree = 0.05;
            this.SpaceReducer = s;
        end
    end
    
    methods(Static)
        
        function [d, r, m] = CreatePlots(dim)
            
            m = models.mathmod2012.MathMODExperiment(dim);
            s = m.SpaceReducer;
            
            d = EstimatorAnalyzer;
            d.EstimatorIterations = [1 2 5];
            d.EstimatorVersions = [1 1 0 0 1 0 0 1 1];
            d.SingleFigures = true;
            
            d.setModel(m);
            r = d.ReducedModel;
            
            for mu1 = [0 1]
                d.start([mu1; 1], 1);

                % 1 = full, 7 = LSLE TD
                md = d.ModelData(end);
                absmax = 30*abs(md.ErrT(7)-md.ErrT(1));
                relmax = 30*abs(md.RelErrT(7)-md.RelErrT(1));

                f = d.Figures{1};
                a = gca(f);
                axis(a,[0 m.T md.MinErr*.9 min(1e4,absmax)]);
                set(legend(a),'Location','NorthEast');
                fi = fullfile(KerMor.App.DataDirectory,'MathMOD',sprintf('MathMOD_mu1_%d_deg%f_', mu1, s.Degree));
                Utils.saveFigure(f,[fi 'errors'],'fig');
                title(a,'');
                Utils.saveFigure(f,[fi 'errors']);

                f = d.Figures{2}; a = gca(f);
                axis(a,[0 m.T md.MinRelErr*.9 min(1,relmax)]);
                set(legend(a),'Location','NorthEast');
                Utils.saveFigure(f,[fi 'relerr'],'fig');
                title(a,'');
                Utils.saveFigure(f,[fi 'relerr']);

                f = d.Figures{3}; a = gca(f);
                Utils.saveFigure(f,[fi 'ctimes'],'fig');
                title(a,'');
                Utils.saveFigure(f,[fi 'ctimes']);
            
            end
            
            PlotParamSweep(r,[1; 1],1,1,-.1:.05:1.1);
            Utils.saveFigure(gcf,'mu1_sweep','png');
            
            d.createStatsTables;
        end
        
    end
end

