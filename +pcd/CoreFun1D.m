classdef CoreFun1D < dscomponents.ACompEvalCoreFun
% The core nonlinear function of the PCD model.
%
% @author Daniel Wirtz @date 2010-03-16
%
% @new{0,6,dw,2011-11-26} Computing the JSparsityPattern
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
%         A;
        
        % The number of nodes
        nodes;
    end
    
    methods
                
        function this = CoreFun1D(dynsys)
            this = this@dscomponents.ACompEvalCoreFun;
            this.sys = dynsys;
            this.MultiArgumentEvaluations = true;
            this.TimeDependent = false;
%             this.CustomJacobian = true;
        end
        
        function copy = clone(this)
            copy = models.pcd.CoreFun1D(this.sys);
            
            % Call superclass method
            copy = clone@dscomponents.ACoreFun(this, copy);
            
            % copy reference!
            %copy.sys = this.sys; % already done in constructor
%             copy.A = this.A;
            copy.nodes = this.nodes;
        end
        
        function newSysDimension(this)
            n = this.sys.Dims(1);
            % Add x_a dependencies
            % 1=x_a, 2=y_a, 3=x_i {, y_i}
            i = [(1:n)'; (1:n)'; (1:n)']; 
            j = (1:3*n)';
            % Add y_a dependencies
            % 1=x_a, 2=y_a {, x_i}, 3=y_i
            i = [i; (n+1:2*n)'; (n+1:2*n)'; (n+1:2*n)'];
            j = [j; (1:2*n)'; (3*n+1:4*n)'];
            % Add x_i dependencies
            % {x_a,} 1=y_a, 2=x_i {, y_i}
            i = [i; (2*n+1:3*n)'; (2*n+1:3*n)']; 
            j = [j; (n+1:3*n)'];
            % Add y_i dependencies
            % 1=x_a, {y_a, x_i}, 2=y_i
            i = [i; (3*n+1:4*n)'; (3*n+1:4*n)']; 
            j = [j; (1:n)'; (3*n+1:4*n)'];
            this.XDim = 4*n;
            this.JSparsityPattern = sparse(i,j,ones(length(i),1),this.XDim,this.XDim);
            this.nodes = n;
        end
        
        function fx = evaluateCoreFun(this, x, ~, mu)
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

                fx(1:m) = mu(1)*xi.*ya - mu(3)*xa + rb;
                fx(m+1:2*m) = mu(2)*yi.*xan - mu(4)*ya;
                fx(2*m+1:3*m) = -mu(1)*xi.*ya - mu(5)*xi + mu(7) - rb;
                fx(3*m+1:end) = -mu(2)*yi.*xan - mu(6)*yi + mu(8);
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

                fx(1:m,:) = bsxfun(@times,xi.*ya,mu(1,:)) - bsxfun(@times,xa,mu(3,:)) + rb;
                fx(m+1:2*m,:) = bsxfun(@times,yi.*xan,mu(2,:)) - bsxfun(@times,ya,mu(4,:));
                fx(2*m+1:3*m,:) = -bsxfun(@times,xi.*ya,mu(1,:)) - bsxfun(@times,xi,mu(5,:)) + bsxfun(@times,ones(size(xi)),mu(7,:)) - rb;
                fx(3*m+1:end,:) = -bsxfun(@times,yi.*xan,mu(2,:)) - bsxfun(@times,yi,mu(6,:)) + bsxfun(@times,ones(size(xi)),mu(8,:));
            end
        end
        
        function fxj = evaluateComponents(this, J, ends, ~, ~, X, ~, mu)
            % The vector embedding results from the fixed ordering of the full 4*m-vector into
            % the components x_a, y_a, x_i, y_i
            m = this.nodes;
            fxj = zeros(length(J),size(X,2));
            
            mu = [repmat(this.sys.ReacCoeff,1,size(mu,2)); mu];
            for idx=1:length(J)
                j = J(idx);
                if idx == 1
                    st = 0;
                else
                    st = ends(idx-1);
                end
                % Select the elements of x that are effectively used in f
                xidx = (st+1):ends(idx);
                x = X(xidx,:);
                
                % X_a
                if j <= m
                    % 1=x_a, 2=y_a, 3=x_i {, y_i}
                    % mu(1)*xi.*ya - mu(3)*xa + rb;
                    fj = mu(1,:).*x(3,:).*x(2,:) - mu(3,:).*x(1,:);
                    % Boundary condition
                    if j == m
                        fj = fj + x(3,:).*mu(9,:)/this.sys.hs;
                    end
                    
                % Y_a
                elseif m < j && j <= 2*m
                    % 1=x_a, 2=y_a {, x_i}, 3=y_i
                    % mu(2)*yi.*xan - mu(4)*ya;
                    fj = mu(2,:).*x(3,:).*x(1,:).^this.sys.n - mu(4,:).*x(2,:);
                    
                % X_i
                elseif 2*m < j && j <= 3*m
                    % {x_a,} 1=y_a, 2=x_i {, y_i}
                    % -mu(1)*xi.*ya - mu(5)*xi + mu(7) - rb;
                    fj = -mu(1,:).*x(2,:).*x(1,:) - mu(5,:).*x(2,:) + mu(7,:);
                    % Boundary condition
                    if j == 3*m
                        fj = fj - x(2,:).*mu(9,:)/this.sys.hs;
                    end
                    
                % Y_i
                else
                    % 1=x_a, {y_a, x_i}, 2=y_i
                    % -mu(2)*yi.*xan - mu(6)*yi + mu(8);
                    fj = -mu(2,:).*x(2,:).*x(1,:).^this.sys.n - mu(6,:).*x(2,:) + mu(8,:);
                end
                fxj(idx,:) = fj;
            end
        end
    end
end

