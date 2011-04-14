classdef CoreFun2D < dscomponents.ACoreFun & ISimConstants
    % The core nonlinear function of the PCD model.
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
        
        function copy = clone(this)
            copy = models.pcd.CoreFun2D(this.sys);
            
            % Call superclass method
            copy = clone@dscomponents.ACoreFun(this, copy);
            
            % copy reference!
            %copy.sys = this.sys; % already done in constructor
            copy.A = this.A;
            copy.dim = this.dim;
            copy.left = this.left;
            copy.right = this.right;
            copy.lower = this.lower;
            copy.upper = this.upper;
        end
        
        function this = CoreFun2D(dynsys)
            this.sys = dynsys;
            this.MultiArgumentEvaluations = true;
        end
        
        function updateSimConstants(this)
            % Create diffusion matrix
            d1 = this.sys.dim1;
            d2 = this.sys.dim2;
            d = d1*d2;
            [this.A,idxmat] = general.MatUtils.laplacemat(this.sys.h,...
                                this.sys.dim1,this.sys.dim2);
            this.upper = idxmat(1:d1:d);
            this.right = idxmat(d-d1+1:d);
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
            
%             if size(x,2) == 1
                % Extract single functions
                xa = x(1:m);
                xan = xa.^n;
                ya = x(m+1:2*m);
                xi = x(2*m+1:3*m);
                yi = x(3*m+1:end);

                % Compile boundary conditions
                rb = zeros(m,1);
%                 rb(this.upper) = -(xi(this.upper)*mu(1))/h;
%                 rb(this.right) = rb(this.right) - (xi(this.right)*mu(2))/h;
%                 rb(this.lower) = rb(this.lower) - (xi(this.lower)*mu(3))/h;
%                 rb(this.left) = rb(this.left) - (xi(this.left)*mu(4))/h;
%                 edges = [this.upper(1) this.upper(end)...
%                          this.lower(1) this.lower(end)];
%                 %rb(edges) = .5*rb(edges);
%                 rb(edges) = 0;

                fx(1:m) = mu(1)*xi.*ya - mu(3)*xa + this.A*xa - rb;
                fx(m+1:2*m) = mu(2)*yi.*xan - mu(4)*ya + this.sys.D2*this.A*ya;
                fx(2*m+1:3*m) = -mu(1)*xi.*ya - mu(5)*xi + mu(7) + this.sys.D3*this.A*xi + rb;
                fx(3*m+1:end) = -mu(2)*yi.*xan - mu(6)*yi + mu(8) + this.sys.D4*this.A*yi;
%             else
%                 % Extract single functions
%                 xa = x(1:m,:);
%                 xan = xa.^n;
%                 ya = x(m+1:2*m,:);
%                 xi = x(2*m+1:3*m,:);
%                 yi = x(3*m+1:end,:);
% 
%                 % Compile boundary conditions
%                 rb = zeros(m,size(x,2));
%                 rb(this.upper,:) = -(bsxfun(@mult,xi(this.upper,:),mu(9,:)))/h;
%                 rb(this.right,:) = rb(this.right,:) - (bsxfun(@mult,xi(this.right,:),mu(10,:)))/h;
%                 rb(this.lower,:) = rb(this.lower,:) - (bsxfun(@mult,xi(this.lower,:),mu(11,:)))/h;
%                 rb(this.left,:) = rb(this.left,:) - (bsxfun(@mult,xi(this.left,:),mu(12,:)))/h;
%                 edges = [this.upper(1) this.upper(end)...
%                          this.lower(1) this.lower(end)];
%                 rb(edges,:) = .5*rb(edges,:);
% 
%                 fx(1:m,:) = bsxfun(@mult,xi.*ya,mu(1,:)) - bsxfun(@mult,xa,mu(3,:)) + this.A*xa - rb;
%                 fx(m+1:2*m,:) = bsxfun(@mult,yi.*xan,mu(2,:)) - bsxfun(@mult,ya,mu(4,:)) + this.sys.D2*this.A*ya;
%                 fx(2*m+1:3*m,:) = -bsxfun(@mult,xi.*ya,mu(1,:)) - bsxfun(@mult,xi,mu(5,:)) + bsxfun(@mult,ones(size(xi)),mu(7,:)) + this.sys.D3*this.A*xi + rb;
%                 fx(3*m+1:end,:) = -bsxfun(@mult,yi.*xan,mu(2,:)) - bsxfun(@mult,yi,mu(6,:)) + bsxfun(@mult,ones(size(xi)),mu(8,:)) + this.sys.D4*this.A*yi;
%             end
            
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

