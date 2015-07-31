classdef BurgersSys_A < models.BaseFirstOrderSystem
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
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties
    end
    
    methods
        function this = BurgersSys_A(model)
            this = this@models.BaseFirstOrderSystem(model);

            % Set core function
            this.f = models.burgers.BurgersF_NoA(this);
            
            this.addParam('v', .5,'Range', [0, 1], 'Desired', 100);
        end
        
        function newDim(this)
            this.updateDimensions;
        end
    end
    
    methods(Access=protected)
        function updateDimensions(this)
            m = this.Model;
            xs  = linspace(m.Omega(1), m.Omega(2), m.Dimension+2)';
            f = @(x)exp(-(15.*(x-.5)).^2);
            x = xs(2:m.Dimension+1);
            x0 = f(x)-f(0);
            this.x0 = dscomponents.ConstInitialValue(x0);
            
            n = m.Dimension;
            dx = (m.Omega(2) - m.Omega(1))/(n+1);
            e = ones(n,1);
            d2 = e/(dx^2);
            
            a = dscomponents.AffLinCoreFun(this);
            a.addMatrix('mu(1)',spdiags([d2 -2*d2  d2], -1:1, n, n));
            this.A = a; 
            
            this.NumStateDofs = n;
            
            this.f.newDim;
            
            updateDimensions@models.BaseFirstOrderSystem(this);
        end
    end
    
end