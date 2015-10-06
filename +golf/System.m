classdef System < models.BaseSecondOrderSystem
% System: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2015-08-06
%
% @new{0,7,dw,2015-08-06} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetAccess=private)
        maxvforhole;
    end
    
    methods
        function this = System(model)
            this = this@models.BaseSecondOrderSystem(model);
            
            % Initial position
            this.addParam('x_0',1);
            this.addParam('y_0',1);
            
            % Initial velocity
            this.addParam('vx_0',1);
            this.addParam('vy_0',0);
            
            % Hole position
            this.addParam('hole_x',0);
            this.addParam('hole_y',0);
            
            % #7: Gravity
            this.addParam('Gravity',9.81);
            % #8: Friction on grass
            this.addParam('Grass friction',0.03);
            
            this.NumStateDofs = 2;
            this.NumDerivativeDofs = 2;
            this.updateDimensions;
            
            x0 = dscomponents.AffineInitialValue;
            x0.addMatrix('mu(1)',[1;0]);
            x0.addMatrix('mu(2)',[0;1]);
            this.x0 = x0;
            
            this.f = models.golf.Force(this);
            
            this.updateSparsityPattern;
        end
        
        function prepareSimulation(this, mu, inputidx)
            prepareSimulation@models.BaseSecondOrderSystem(this, mu, inputidx);
            m = this.Model;
            
            maxv = sqrt((1/(2*m.rad_ball)) * mu(7) * (2*m.rad_hole - m.rad_ball)^2 ); % m/s
            this.maxvforhole = maxv;
            
            % Evtl. den Zeitschritt verkleinern, so dass der Ball nicht pro Schritt
            % die Lochweite überspringen kann ODER
            % Mindestens eine Auflösung von 24 Bildern pro sekunde, damit beim
            % Aviexport ein flüssiges Video entsteht!
            dt = m.dt;
            if (dt > 2*m.rad_hole/maxv) || (dt > 1/24)
                this.MaxTimestep = min(0.9*2*m.rad_hole/maxv, 1/24);
            else
                this.MaxTimestep = dt;
            end
        end
        
        function dv = getDerivativeDirichletValues(this)
            dv = [];
        end
    end
    
    methods(Access=protected)
        function val = getDerivativeInitialValues(~, mu)
            val = mu(3:4);
        end
    end
    
end