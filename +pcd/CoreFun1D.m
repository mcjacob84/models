classdef CoreFun1D < dscomponents.ACoreFun
% The core nonlinear function of the PCD model.
%
% @author Daniel Wirtz @date 2010-03-16
%
% @change{0,5,dw,2011-11-02} Augmenting the mu parameters by the base system's
% models.pcd.BasePCDSystem.ReacCoeff vector. This removes the reaction coefficients from
% the system as true parameters but allows to quickly revert the process if needed.
%
% @change{0,3,dw,2011-04-12} Added MultiArgumentEvaluation capabilities to this class.
%
% @new{0,3,dw,2011-03-16} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(Access=private)
        % The assoc. dynamical system
        % (need some values from that)
        sys;
    
        % The diffusion matrix
        A;
        
        % The number of nodes
        nodes;
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
            copy.nodes = this.nodes;
        end
        
        function newSysDimension(this)
            % Create diffusion matrix
            d = this.sys.Dims(1);
            e = ones(d,1);
            this.A = spdiags([e -2*e e],-1:1,d,d)/this.sys.hs^2;
            this.A(1,2) = 2/this.sys.hs^2;
            this.A(end,end-1) = 2/this.sys.hs^2;
            this.nodes = this.sys.Dims(1);
        end
        
        function fx = evaluateCoreFun(this, x, t, mu)%#ok
            % Allocate result vector
            fx = zeros(size(x));
            
            m = this.nodes;
            
            % Non-Multidimensional case
            if size(x,2) == 1
                mu = [this.sys.ReacCoeff; mu];
                
                % Extract single functions
                xa = x(1:m);
                xan = xa.^this.sys.n;
                ya = x(m+1:2*m);
                xi = x(2*m+1:3*m);
                yi = x(3*m+1:end);

                % Boundary conditions
                rb = zeros(m,1);
                rb(end) = xi(end)*mu(9)/this.sys.hs; %max((500-t)/500,0)*

                fx(1:m) = mu(1)*xi.*ya - mu(3)*xa + this.A*xa + rb;
                fx(m+1:2*m) = mu(2)*yi.*xan - mu(4)*ya + this.sys.Diff(1)*this.A*ya;
                fx(2*m+1:3*m) = -mu(1)*xi.*ya - mu(5)*xi + mu(7) + this.sys.Diff(2)*this.A*xi - rb;
                fx(3*m+1:end) = -mu(2)*yi.*xan - mu(6)*yi + mu(8) + this.sys.Diff(3)*this.A*yi;
            else
                mu = [repmat(this.sys.ReacCoeff,1,size(mu,2)); mu];
                % Extract single functions
                xa = x(1:m,:);
                xan = xa.^this.sys.n;
                ya = x(m+1:2*m,:);
                xi = x(2*m+1:3*m,:);
                yi = x(3*m+1:end,:);

                % Boundary conditions
                rb = zeros(m,size(xi,2));
                rb(end,:) = (xi(end,:).*mu(9,:))/this.sys.hs;

                fx(1:m,:) = bsxfun(@times,xi.*ya,mu(1,:)) - bsxfun(@times,xa,mu(3,:)) + this.A*xa + rb;
                fx(m+1:2*m,:) = bsxfun(@times,yi.*xan,mu(2,:)) - bsxfun(@times,ya,mu(4,:)) + this.sys.Diff(1)*this.A*ya;
                fx(2*m+1:3*m,:) = -bsxfun(@times,xi.*ya,mu(1,:)) - bsxfun(@times,xi,mu(5,:)) + bsxfun(@times,ones(size(xi)),mu(7,:)) + this.sys.Diff(2)*this.A*xi - rb;
                fx(3*m+1:end,:) = -bsxfun(@times,yi.*xan,mu(2,:)) - bsxfun(@times,yi,mu(6,:)) + bsxfun(@times,ones(size(xi)),mu(8,:)) + this.sys.Diff(3)*this.A*yi;
            end
        end
    end
end

