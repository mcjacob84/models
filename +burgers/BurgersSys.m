classdef BurgersSys < models.BaseDynSystem
% BurgersSys: 
%
%
%
% @author Daniel Wirtz @date 2012-04-24
%
% @new{0,6,dw,2012-04-24} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
    end
    
    methods
        function this = BurgersSys(model)
            this = this@models.BaseDynSystem(model);
            
            % Set core function
            this.f = models.burgers.BurgersF(this);
            
            this.addParam('v', [0, 1], 100);
        end
        
        function newDim(this)
            m = this.Model;
            xs  = linspace(m.Omega(1), m.Omega(2), m.Dimension+2)';
            f = inline('exp(-(15.*(x-.5)).^2)','x');
            x = xs(2:m.Dimension+1);
            x0 = f(x)-f(0);
            this.x0 = dscomponents.ConstInitialValue(x0);
            this.f.newDim;
        end
    end
    
end