classdef Belly < models.muscle.AMuscleConfig
    
    properties(Access=private)
        NumParts;
    end
    
    methods
        function this = Belly(varargin)
            this = this@models.muscle.AMuscleConfig(varargin{:});
            this.init;
            
            %% Muscle fibre weights
            types = [0 .2 .4 .6 .8 1];
            ftw = zeros(this.FEM.GaussPointsPerElem,length(types),this.FEM.Geometry.NumElements);
            % Test: Use only slow-twitch muscles
            ftw(:,1,:) = .4;
            ftw(:,2,:) = .05;
            ftw(:,3,:) = .05;
            ftw(:,4,:) = .1;
            ftw(:,5,:) = .2;
            ftw(:,6,:) = .2;
            this.FibreTypeWeights = ftw;
            p = models.motorunit.Pool;
            p.FibreTypes = types;
            this.Pool = p;
        end
        
        function configureModel(this, m)
            configureModel@models.muscle.AMuscleConfig(this, m);
            m.T = 150;
            m.dt = .1;
            m.DefaultMu(4) = 6;
            m.Plotter.DefaultArgs = {'Pool',true};
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(this)
            np = 4;
            geo = fem.geometry.Belly(np,10,'Radius',1,'InnerRadius',.5,'Gamma',2);
            this.NumParts = np;
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.FEM.Geometry;
            for k = geo.NumElements-3:geo.NumElements
                displ_dir(:,geo.Elements(k,geo.MasterFaces(4,:))) = true;
            end
            for k = 1:4
                displ_dir(:,geo.Elements(k,geo.MasterFaces(3,:))) = true;
            end
        end
        
        function anull = seta0(~, anull)
            anull(2,:,:) = 1;
        end
    end
    
    methods(Static)
        function res = test_BellyGeometryGeneration
            g = fem.geometry.Belly(4,35,'InnerRadius',.2,'Gamma',7);
            g = fem.geometry.Belly(4,35,'InnerRadius',[.2 .6],'Gamma',5);
            g = fem.geometry.Belly(4,35,'InnerRadius',.2,'Gamma',[10 20]);
            g = fem.geometry.Belly(4,35,'InnerRadius',[.2 .6],'Gamma',[10 20]);
            g = fem.geometry.Belly(4,35,'Radius',[4 2],'InnerRadius',.5,'Gamma',[10 20]);
            g = fem.geometry.Belly(4,35,'Radius',@(x)[sqrt(abs(x)); 1./(x-34).^2],'InnerRadius',.2);
            g.plot;
            g = fem.geometry.Belly(4,35,'Radius',@(x)sqrt(abs(x)));
            g.plot;
%             g = fem.geometry.Belly(4,35,'Radius',[4 2],'InnerRadius',.5,...
%                 'Gamma',[10 20],'MinZ',-1.5);
%             g.plot;
            g = fem.geometry.Belly(4,35,'Radius',@(x)sqrt(abs(x)),'Layers',[.4 .7 1]);
            g.plot;
            res = 1;
        end
        
        function test_BellyModel
            m = models.muscle.Model(models.muscle.examples.Belly);
            m.simulateAndPlot;
        end
    end
    
end
