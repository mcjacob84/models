classdef CubePull < models.muscle.AMuscleConfig
    
    methods
        function this = CubePull(varargin)
            % Single cube with same config as reference element
            this = this@models.muscle.AMuscleConfig(varargin{:});
            this.init;
        end
        
        function configureModel(this, m)
            configureModel@models.muscle.AMuscleConfig(this, m);
            m.T = 150;
            m.dt = .5;
            m.DefaultMu(1) = .5;
            m.DefaultMu(3) = .1;
            m.DefaultInput = 1;
            
            % MOR pre-setup.
            % If you assign a new SpaceReducer instance, dont forget to set
            % the TargetDimensions property accordingly (is now
            % automatically set within the configureModel base routine)
            %s = spacereduction.PODGreedy;
            %s.Eps = 1e-9;
            %s.MaxSubspaceSize = 500;
            s = spacereduction.PODReducer;
            s.IncludeInitialSpace = true;
            s.IncludeBSpan = true;
            s.Mode = 'abs';
            m.SpaceReducer = s;
        end
        
        function configureModelFinal(this)
            % Set desired reduction to full state space dimension by
            % default
            m = this.Model;
            m.SpaceReducer.Value = m.System.NumDerivativeDofs;
            configureModelFinal@models.muscle.AMuscleConfig(this);
        end
        
        function P = getBoundaryPressure(this, elemidx, faceidx)
            % Determines the neumann forces on the boundary.
            %
            % The unit for the applied quantities is megaPascal [MPa]
            %
            % In the default implementation there are no force boundary
            % conditions.
            P = [];
            if this.Options.GeoNr == 2
                if any(elemidx == [1 2 61 62]) && faceidx == 3
                    P = 1;
                end
            else
                if elemidx == 1 && faceidx == 3
                    P = 1;
                end
            end
        end
        
        function u = getInputs(this)
            u = {this.getAlphaRamp(this.Model.T/3,1)};
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(this)
            if this.Options.GeoNr == 2
                geo = fem.geometry.RegularHex27Grid([0 2.5 5],-10:20,[0 2.5 5]);
            else
                geo = fem.geometry.RegularHex27Grid([0 5],-1:10,[0 5]);
            end
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.FEM.Geometry;
            % Fix all on left and only the y,z directions of the back right
            % nodes
            if this.Options.GeoNr == 2
                displ_dir(2,geo.Elements([59 60 119 120],geo.MasterFaces(4,:))) = true;
            else
                displ_dir(2,geo.Elements(end,geo.MasterFaces(4,:))) = true;
            end
            % Fix middle point (works for both geo versions)
            displ_dir(:,geo.Elements(end,geo.MasterFaces(4,1))) = true;
        end
        
        function anull = seta0(~, anull)
           % Direction is xz
           %anull([1 3],:,:) = 1;
        end
    end
    
    methods(Static)
        function test_CubePull
            m = models.muscle.Model(models.muscle.examples.CubePull);
            mu = m.getRandomParam;
            mu(3) = .1;
            m.simulateAndPlot(true,mu,1);
        end
    end
end

