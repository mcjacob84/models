classdef CoreFun1D < dscomponents.ACoreFun & ISimConstants
    %CoreFun The core nonlinear function of the PCD model.
    %
    % @author Daniel Wirtz @date 16.03.2010
    
    properties(Access=private)
        % The assoc. dynamical system
        % (need some values from that)
        sys;
        
        A;
    end
    
    methods
        
        function res = project(this, V)%#ok
            error('Projection of pure model core function not supported in this model.');
        end
        
        function this = CoreFun1D(dynsys)
            this.sys = dynsys;
        end
        
        function updateSimConstants(this)
            % Create diffusion matrix
            d = this.sys.dim;
            e = ones(d,1);
            this.A = spdiags([e -2*e e],-1:1,d,d);
        end
        
        function fx = evaluateCoreFun(this, x, t, mu)%#ok
            % Allocate result vector
            fx = zeros(size(x));
            
            m = this.sys.dim;
            
            % Extract single functions
            xa = x(1:m);
            xan = xa.^this.sys.n;
            ya = x(m+1:2*m);
            xi = x(2*m+1:3*m);
            yi = x(3*m+1:end);
            
            % Boundary conditions
            rb = zeros(m,1);
            rb(end) = - (xi(end)*mu(1))/this.sys.h;
            
            fx(1:m) = mu(2)*this.sys.lam1*xi.*ya - mu(4)*xa + this.A*xa - rb;
            fx(m+1:2*m) = mu(3)*this.sys.lam2*yi.*xan - mu(5)*ya + this.sys.D*this.A*ya;
            fx(2*m+1:3*m) = -mu(2)*xi.*ya - mu(4)*xi + mu(6) + this.A*xi + rb;
            fx(3*m+1:end) = -mu(3)*yi.*xan - mu(5)*yi + mu(7) + this.sys.D*this.A*yi;
        end
        
        function c = getGlobalLipschitz(this, mu, inputidx)
            c = 1;
        end
    end
end

