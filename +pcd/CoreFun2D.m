classdef CoreFun2D < dscomponents.ACoreFun & ISimConstants
    %CoreFun The core nonlinear function of the PCD model.
    %
    % @author Daniel Wirtz @date 16.03.2010
    
    properties(Access=private)
        % The assoc. dynamical system
        % (need some values from that)
        sys;
        
        upper;
        lower;
        right;
        left;
        A;
        dim;
    end
    
    methods
        
        function res = project(this, V)%#ok
            error('Projection of pure model core function not supported in this model.');
        end
        
        function this = CoreFun2D(dynsys)
            this.sys = dynsys;
        end
        
        function updateSimConstants(this)
            % Create diffusion matrix
            d1 = this.sys.dim1;
            d2 = this.sys.dim2;
            d = d1*d2;
            [this.A,idxmat] = general.MatUtils.laplacemat(this.sys.h,...
                                this.sys.dim1,this.sys.dim2);
            this.upper = idxmat(1:d1:d);
            this.right = idxmat(d-d1:d);
            this.lower = idxmat(d1:d1:d);
            this.left = idxmat(1:d1);
            
            this.dim = d;
        end
        
        function fx = evaluateCoreFun(this, x, t, mu)%#ok
            % Allocate result vector
            fx = zeros(size(x));
            
            m = this.dim;
            n = this.sys.n;
            h = this.sys.h;
            
            % Extract single functions
            xa = x(1:m);
            xan = xa.^n;
            ya = x(m+1:2*m);
            xi = x(2*m+1:3*m);
            yi = x(3*m+1:end);
            
            % Compile boundary conditions
            rb = zeros(m,1);
            rb(this.upper) = -(xi(this.upper)*mu(1))/h;
            rb(this.right) = rb(this.right) - (xi(this.right)*mu(2))/h;
            rb(this.lower) = rb(this.lower) - (xi(this.lower)*mu(3))/h;
            rb(this.left) = rb(this.left) - (xi(this.left)*mu(4))/h;
            edges = [this.upper(1) this.upper(end)...
                     this.lower(1) this.lower(end)];
            rb(edges) = .5*rb(edges);
            
            % Handle xa function
            %         fx(1:m) = A*xa;
            %         fx(m+1:2*m) = A*ya;
            %         fx(2*m+1:3*m) = A*xi;
            %         fx(3*m+1:end) = A*yi;
            fx(1:m) = mu(5)*this.sys.lam1*xi.*ya - mu(7)*xa + this.A*xa - rb;
            fx(m+1:2*m) = mu(6)*this.sys.lam2*yi.*xan - mu(8)*ya + this.sys.D*this.A*ya;
            fx(2*m+1:3*m) = -mu(5)*xi.*ya - mu(7)*xi + mu(9) + this.A*xi + rb;
            fx(3*m+1:end) = -mu(6)*yi.*xan - mu(8)*yi + mu(10) + this.sys.D*this.A*yi;
            
            %         figure(1);
            %         plo(x,1);
            %         title('casp-8');
            %         plo(x,2);
            %         title('casp-3');
            %         plo(x,3);
            %         title('procasp-8');
            %         plo(x,4);
            %         title('procasp-3');
            %         pause;
        end
    end
    
%         function plo(fx,idx)
%         subplot(2,2,idx);
%         x = fx((idx-1)*m+1:idx*m);
%         surf(reshape(x,d1,d2));
%     end
end

