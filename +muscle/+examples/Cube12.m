classdef Cube12 < models.muscle.AMuscleConfig
    
    methods
        function this = Cube12(varargin)
            this = this@models.muscle.AMuscleConfig(varargin{:});
            this.init;
        end
        
        function configureModel(this, m)
            configureModel@models.muscle.AMuscleConfig(this, m);
            m.T = 40;
            m.dt = .05;
            % Activate over 20ms
            m.DefaultMu(2) = 20;
            m.DefaultMu(13) = .1; % [MPa]
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(~)
            geo = fem.geometry.RegularHex8Grid(-1:1,-1:2,-1:1);
            geo = geo.toCube20Node;
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            geo = this.FEM.Geometry;
            for k = [5 6 11 12]
                displ_dir(:,geo.Elements(k,geo.MasterFaces(4,:))) = true;
            end
        end

        function anull = seta0(~, anull)
            % Direction is xz
            anull([1 3],:,:) = 1;
        end
    end
    
    methods(Static)
        function test_Cube12
            m = models.muscle.Model(models.muscle.examples.Cube12);
            m.simulateAndPlot;
        end
    end
    
end

