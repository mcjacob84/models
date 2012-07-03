classdef CoreFun2D < dscomponents.ACompEvalCoreFun
% The core nonlinear function of the PCD model.
%
% @author Daniel Wirtz @date 16.03.2010
%
% @change{0,5,dw,2011-11-02} Augmenting the mu parameters by the base system's
% models.pcd.BasePCDSystem.ReacCoeff vector. This removes the reaction coefficients from
% the system as true parameters but allows to quickly revert the process if needed.
%
% @change{0,5,dw,2011-10-17} Removed the ISimConstants class and
% unified the structure of the 1D-3D pcd models.
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
        
        nodes;
        
        hlp;
        
        idxmat;
    end
    
    methods
        
        function this = CoreFun2D(dynsys)
            this = this@dscomponents.ACompEvalCoreFun;
            this.sys = dynsys;
            this.MultiArgumentEvaluations = true;
            this.TimeDependent = false;
            this.hlp.n = this.sys.n;
        end
        
        function copy = clone(this)
            % sys already copied in constructor (see below)
            copy = models.pcd.CoreFun2D(this.sys);
            
            % Call superclass method
            copy = clone@dscomponents.ACoreFun(this, copy);
            
            copy.nodes = this.nodes;
        end
        
        function newSysDimension(this)
            % Create diffusion matrix
            s = this.sys;            
            n = prod(s.Dims);
            this.nodes = n;
            this.hlp.d1 = s.Dims(1);
            this.hlp.d2 = s.Dims(2);
            
            this.idxmat = zeros(s.Dims(1),s.Dims(2));
            this.idxmat(:) = 1:this.nodes;
            
            % Can use the unscaled h and original Omega for computation of
            % ranges & distances from middle
            this.hlp.hs = s.hs;
            a = s.Omega;
            this.hlp.xr = a(1,2)-a(1,1); % xrange
            this.hlp.yr = a(2,2)-a(2,1); % yrange
            this.hlp.xd = abs(((1:this.hlp.d1)-1)*this.sys.h-.5*this.hlp.xr); % x distances from middle
            this.hlp.yd = abs(((1:this.hlp.d2)-1)*this.sys.h-.5*this.hlp.yr); % y distances from middle
            
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
        end
        
        function fx = evaluateCoreFun(this, x, ~, mu)
            % Allocate result vector
            fx = zeros(size(x));
            
            m = this.nodes;
            
            % Compile boundary conditions
            %to = this.sys.Model.tau*t;
            %ud = (to < 10)*1 + (to >= 10)*max(0,(2-to/10));
            ud = 1;
            
            if size(x,2) == 1
                % Uncomment if reaction coeffs become real params again
                %mu = [s.ReacCoeff; mu]';
                mu = [this.sys.ReacCoeff; mu([1 1 1 1 2 2 2 2])];
                
                % Extract single functions
                xa = x(1:m);
                xan = xa.^this.hlp.n;
                ya = x(m+1:2*m);
                xi = x(2*m+1:3*m);
                yi = x(3*m+1:end);
    
                rb = zeros(m,1);
                
                %% Top & Bottom
                % bottom
                pos = this.hlp.xd <= this.hlp.xr*mu(10)/2;
                idx = this.idxmat(pos,1); 
                rb(idx) = (xi(idx)*mu(14)*ud);
                % top
                pos = this.hlp.xd <= this.hlp.xr*mu(9)/2;
                idx = this.idxmat(pos,end);
                rb(idx) = rb(idx) + (xi(idx)*mu(13)*ud);
                
                %% Left & Right
                % right
                pos = this.hlp.yd <= this.hlp.yr*mu(12)/2;
                idx = this.idxmat(end,pos);
                rb(idx) = rb(idx) + (xi(idx)*mu(16)*ud);
                % left
                pos = this.hlp.yd <= this.hlp.yr*mu(11)/2;
                idx = this.idxmat(1,pos);
                rb(idx) = rb(idx) + (xi(idx)*mu(15)*ud);
                
                % x_a
                fx(1:m) = mu(1)*xi.*ya - mu(3)*xa + rb/this.hlp.hs;
                % y_a
                fx(m+1:2*m) = mu(2)*yi.*xan - mu(4)*ya;
                % x_i
                fx(2*m+1:3*m) = -mu(1)*xi.*ya - mu(5)*xi + mu(7) - rb/this.hlp.hs;
                % y_i
                fx(3*m+1:end) = -mu(2)*yi.*xan - mu(6)*yi + mu(8);
            else
                % Uncomment if reaction coeffs become real params again
                %mu = [s.ReacCoeff; mu]';
                mu = [repmat(this.sys.ReacCoeff,1,size(mu,2)); mu([1 1 1 1 2 2 2 2],:)];
                
                % Extract single functions
                xa = x(1:m,:);
                xan = xa.^this.hlp.n;
                ya = x(m+1:2*m,:);
                xi = x(2*m+1:3*m,:);
                yi = x(3*m+1:end,:);

                % Compile boundary conditions
                nd = size(x,2);
                rb = zeros(m,nd);
                
                %% Top & Bottom
                xd = repmat(this.hlp.xd',1,nd); % y distances
                % bottom
                pos = bsxfun(@lt,xd,this.hlp.xr*mu(10,:)/2);
                idx = this.idxmat(:,1);
                rb(idx,:) = pos .* (bsxfun(@mult,xi(idx,:),mu(14,:)*ud));
                % top
                pos = bsxfun(@lt,xd,this.hlp.xr*mu(9,:)/2);
                idx = this.idxmat(:,end);
                rb(idx,:) = rb(idx,:) + pos .* (bsxfun(@mult,xi(idx,:),mu(13,:)*ud));
                
                %% Left & Right
                yd = repmat(this.hlp.yd',1,nd);
                % right
                pos = bsxfun(@lt,yd,this.hlp.yr*mu(12,:)/2);
                idx = this.idxmat(end,:);
                rb(idx,:) = rb(idx,:) + pos .* (bsxfun(@mult,xi(idx,:),mu(16,:)*ud));
                % left
                pos = bsxfun(@lt,yd,this.hlp.yr*mu(11,:)/2);
                idx = this.idxmat(1,:);
                rb(idx,:) = rb(idx,:) + pos .* (bsxfun(@mult,xi(idx,:),mu(15,:)*ud));
                
                fx(1:m,:) = bsxfun(@mult,xi.*ya,mu(1,:)) - bsxfun(@mult,xa,mu(3,:)) + rb/this.hlp.hs;
                fx(m+1:2*m,:) = bsxfun(@mult,yi.*xan,mu(2,:)) - bsxfun(@mult,ya,mu(4,:));
                fx(2*m+1:3*m,:) = -bsxfun(@mult,xi.*ya,mu(1,:)) - bsxfun(@mult,xi,mu(5,:)) + bsxfun(@mult,ones(size(xi)),mu(7,:)) - rb/this.hlp.hs;
                fx(3*m+1:end,:) = -bsxfun(@mult,yi.*xan,mu(2,:)) - bsxfun(@mult,yi,mu(6,:)) + bsxfun(@mult,ones(size(xi)),mu(8,:));
            end
        end
        
        function fxj = evaluateComponents(this, J, ends, ~, ~, X, ~, mu)
            % The vector embedding results from the fixed ordering of the full 4*m-vector into
            % the components x_a, y_a, x_i, y_i
            m = this.nodes;
            nd = size(X,2);
            fxj = zeros(length(J),nd);
            
            %mu = [repmat(this.sys.ReacCoeff,1,size(mu,2)); mu];
            mu = [repmat(this.sys.ReacCoeff,1,size(mu,2)); mu([1 1 1 1 2 2 2 2],:)];
            if nd > 1
                xd = repmat(this.hlp.xd',1,nd);
                yd = repmat(this.hlp.yd',1,nd);
                bottom = bsxfun(@lt,xd,this.hlp.xr*mu(10,:)/2);
                top = bsxfun(@lt,xd,this.hlp.xr*mu(9,:)/2);
                right = bsxfun(@lt,yd,this.hlp.yr*mu(12,:)/2);
                left = bsxfun(@lt,yd,this.hlp.yr*mu(11,:)/2);
            else
                bottom = this.hlp.xd <= this.hlp.xr*mu(10)/2;
                top = this.hlp.xd <= this.hlp.xr*mu(9)/2;
                right = this.hlp.yd <= this.hlp.yr*mu(12)/2;
                left = this.hlp.yd <= this.hlp.yr*mu(11)/2;
            end
            % Get matrix indices
            J2 = mod(J,m);
            J2(J2==0) = m;
            [row, col] = ind2sub([this.hlp.d1 this.hlp.d2],J2);
            
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
                    
                    % Boundary conditions
                    % Bottom
                    if col(idx) == 1
                        fj = fj + bottom(row(idx),:) .* (x(3,:).*mu(14,:)/this.hlp.hs);
                    end
                    % Top
                    if col(idx) == this.hlp.d2
                        fj = fj + top(row(idx),:) .* (x(3,:).*mu(13,:)/this.hlp.hs);
                    end
                    % Right
                    if row(idx) == this.hlp.d1
                        fj = fj + right(col(idx),:) .* (x(3,:).*mu(16,:)/this.hlp.hs);
                    end
                    % Left
                    if row(idx) == 1
                        fj = fj + left(col(idx),:) .* (x(3,:).*mu(15,:)/this.hlp.hs);
                    end
                    
                % Y_a
                elseif m < j && j <= 2*m
                    % 1=x_a, 2=y_a {, x_i}, 3=y_i
                    % mu(2)*yi.*xan - mu(4)*ya;
                    fj = mu(2,:).*x(3,:).*x(1,:).^this.hlp.n - mu(4,:).*x(2,:);
                    
                % X_i
                elseif 2*m < j && j <= 3*m
                    % {x_a,} 1=y_a, 2=x_i {, y_i}
                    % -mu(1)*xi.*ya - mu(5)*xi + mu(7) - rb;
                    fj = -mu(1,:).*x(2,:).*x(1,:) - mu(5,:).*x(2,:) + mu(7,:);
                    % Bottom
                    if col(idx) == 1
                        fj = fj - bottom(row(idx),:) .* (x(2,:).*mu(14,:)/this.hlp.hs);
                    end
                    % Top
                    if col(idx) == this.hlp.d2
                        fj = fj - top(row(idx),:) .* (x(2,:).*mu(13,:)/this.hlp.hs);
                    end
                    % Right
                    if row(idx) == this.hlp.d1
                        fj = fj - right(col(idx),:) .* (x(2,:).*mu(16,:)/this.hlp.hs);
                    end
                    % Left
                    if row(idx) == 1
                        fj = fj - left(col(idx),:) .* (x(2,:).*mu(15,:)/this.hlp.hs);
                    end
                    
                % Y_i
                else
                    % 1=x_a, {y_a, x_i}, 2=y_i
                    % -mu(2)*yi.*xan - mu(6)*yi + mu(8);
                    fj = -mu(2,:).*x(2,:).*x(1,:).^this.hlp.n - mu(6,:).*x(2,:) + mu(8,:);
                end
                fxj(idx,:) = fj;
            end
        end
    end
end

