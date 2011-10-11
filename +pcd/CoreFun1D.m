classdef CoreFun1D < dscomponents.ACoreFun & ISimConstants
    % The core nonlinear function of the PCD model.
    %
    % @author Daniel Wirtz @date 16.03.2010
    %
    % @change{0,3,dw,2011-04-12} Added MultiArgumentEvaluation capabilities to this class.
    
    properties(Access=private)
        % The assoc. dynamical system
        % (need some values from that)
        sys;
    end
    
    properties(Transient)
        % The diffusion matrix
        A;
    end
    
    methods
                
        function this = CoreFun1D(dynsys)
            this.sys = dynsys;
            this.MultiArgumentEvaluations = true;
        end
        
        function copy = clone(this)
            copy = models.pcd.CoreFun1D(this.sys);
            
            % Call superclass method
            copy = clone@dscomponents.ACoreFun(this, copy);
            
            % copy reference!
            %copy.sys = this.sys; % already done in constructor
            copy.A = this.A;
        end
        
        function prepareConstants(this)
            % Create diffusion matrix
            d = this.sys.dim;
            e = ones(d,1);
            this.A = spdiags([e -2*e e],-1:1,d,d)/this.sys.hs^2;
            this.A(1,2) = 2/this.sys.hs^2;
            this.A(end,end-1) = 2/this.sys.hs^2;
        end
        
        function fx = evaluateCoreFun(this, x, t, mu)%#ok
            % Allocate result vector
            fx = zeros(size(x));
            
            m = this.sys.dim;
            
            % Non-Multidimensional case
            if size(x,2) == 1
                % Extract single functions
                xa = x(1:m);
                xan = xa.^this.sys.n;
                ya = x(m+1:2*m);
                xi = x(2*m+1:3*m);
                yi = x(3*m+1:end);

                % Boundary conditions
                rb = zeros(m,1);
%                 if xi(end) ~= 0
%                    keyboard;
%                 end
                rb(end) = - (xi(end)*mu(9))/this.sys.hs;
                %rb(end) = -1/this.sys.hs;
                %rb(end) = 1;

                fx(1:m) = mu(1)*xi.*ya - mu(3)*xa + this.A*xa - rb;
                fx(m+1:2*m) = mu(2)*yi.*xan - mu(4)*ya + this.sys.D2*this.A*ya;
                fx(2*m+1:3*m) = -mu(1)*xi.*ya - mu(5)*xi + mu(7) + this.sys.D3*this.A*xi + rb;
                fx(3*m+1:end) = -mu(2)*yi.*xan - mu(6)*yi + mu(8) + this.sys.D4*this.A*yi;
                %plot(fx); pause;
            else
                % Extract single functions
                xa = x(1:m,:);
                xan = xa.^this.sys.n;
                ya = x(m+1:2*m,:);
                xi = x(2*m+1:3*m,:);
                yi = x(3*m+1:end,:);

                % Boundary conditions
                rb = zeros(m,size(xi,2));
                rb(end,:) = - (xi(end,:).*mu(9,:))/this.sys.hs;

                fx(1:m,:) = bsxfun(@times,xi.*ya,mu(1,:)) - bsxfun(@times,xa,mu(3,:)) + this.A*xa - rb;
                fx(m+1:2*m,:) = bsxfun(@times,yi.*xan,mu(2,:)) - bsxfun(@times,ya,mu(4,:)) + this.sys.D2*this.A*ya;
                fx(2*m+1:3*m,:) = -bsxfun(@times,xi.*ya,mu(1,:)) - bsxfun(@times,xi,mu(5,:)) + bsxfun(@times,ones(size(xi)),mu(7,:)) + this.sys.D3*this.A*xi + rb;
                fx(3*m+1:end,:) = -bsxfun(@times,yi.*xan,mu(2,:)) - bsxfun(@times,yi,mu(6,:)) + bsxfun(@times,ones(size(xi)),mu(8,:)) + this.sys.D4*this.A*yi;
            end
        end
        
        function c = getGlobalLipschitz(this, mu, inputidx)%#ok
            c = 1;
        end
    end
end

