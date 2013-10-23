classdef InhibitCoreFun2D < models.pcdi.CoreFun2D
% The core nonlinear function of the PCD model.
%
% @author Daniel Wirtz @date 16.03.2010
%
% @new{0,6,dw,2012-07-16} Added direct jacobian evaluation function.
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
    
    methods
        
        function this = InhibitCoreFun2D(dynsys)
            this = this@models.pcdi.CoreFun2D(dynsys);
        end
           
        function copy = clone(this)
            % System already copied in constructor (see below)
            copy = models.pcdi.InhibitCoreFun2D(this.System);
            
            % Call superclass method
            copy = clone@models.pcdi.CoreFun2D(this, copy);
        end
        
        function newSysDimension(this)
            % Call basic stuff
            newSysDimension@models.pcdi.CoreFun2D(this);
            
            % Re-build sparsity pattern
            n = this.nodes;
            
            % 1=x_a, 2=y_a, 3=x_i, 4=y_i, 5=iap, 6=bar, 7=yb, 8=xb
            % Add x_a dependencies
            % rc(2)*xi.*ya - rc(5)*xa -rc(11)*xa.*bar + rc(18)*xb
            i = repmat(pos(1),1,5); 
            j = pos([1:3 6 8]);
            % Add y_a dependencies
            % rc(1)*yi.*xa - rc(6)*ya - rc(3)*ya.*iap + rc(14)*yb;
            i = [i repmat(pos(2),1,5)];
            j = [j pos([1:2 4 5 7])];
            % Add x_i dependencies
            % -rc(2)*xi.*ya - rc(9)*xi + rc(16) - rb/this.hlp.hs;
            i = [i repmat(pos(3),1,2)]; 
            j = [j pos([2 3])];
            % Add y_i dependencies
            % -rc(1)*yi.*xa - rc(10)*yi + rc(17);
            i = [i repmat(pos(4),1,2)]; 
            j = [j pos([1 4])];
            
            % 1=x_a, 2=y_a, 3=x_i, 4=y_i, 5=iap, 6=bar, 7=yb, 8=xb
            % Add iap dependencies
            % -(rc(3)+rc(4))*ya.*iap - rc(8)*iap + rc(15) + rc(14)*yab;
            i = [i repmat(pos(5),1,3)]; 
            j = [j pos([2 5 7])];
            % Add bar dependencies
            % -rc(11)*xa.*bar + rc(18)*xab - rc(19)*bar + rc(19);
            i = [i repmat(pos(6),1,3)]; 
            j = [j pos([1 6 8])];
            % Add bar dependencies
            % rc(3)*ya.*iap -(rc(14)+rc(7))*yab;
            i = [i repmat(pos(7),1,3)]; 
            j = [j pos([2 5 7])];
            % Add bar dependencies
            % rc(11)*xa.*bar-(rc(18)+rc(13))*xab;
            i = [i repmat(pos(8),1,3)]; 
            j = [j pos([1 6 8])];
            this.xDim = 8*n;
            this.fDim = 8*n;
            this.JSparsityPattern = sparse(i,j,ones(length(i),1),this.fDim,this.xDim);
            
            function idx = pos(nr)
                idx = [];
                for k=1:length(nr)
                    idx = [idx (nr(k)-1)*n+1:nr(k)*n];%#ok
                end
            end
        end
        
        function fx = evaluate(this, x, t, mu)
            % If this has been projected, restore full size and compute values.
            if ~isempty(this.V)
                x = this.V*x;
            end
            
            % Allocate result vector
            fx = zeros(size(x));

            m = this.nodes;
            
            % Compute activation fun values (it's set up in scaled times!)
            ud = this.activationFun(t,mu);
            rc = this.System.ReacCoeff;
            if size(x,2) == 1
               
                % Extract single functions
                xa = x(1:m);
                xan = xa.^mu(4);
                ya = x(m+1:2*m);
                xi = x(2*m+1:3*m);
                yi = x(3*m+1:4*m);
                iap = x(4*m+1:5*m);
                bar = x(5*m+1:6*m);
                yab = x(6*m+1:7*m);
                xab = x(7*m+1:end);
    
                rb = zeros(m,1);
                
                %% Top & Bottom
                % bottom
                pos = this.hlp.xd <= this.hlp.xr*mu(1)/2;
                idx = this.idxmat(pos,1); 
                rb(idx) = (xi(idx)*mu(2)*ud);
                % top
                pos = this.hlp.xd <= this.hlp.xr*mu(1)/2;
                idx = this.idxmat(pos,end);
                rb(idx) = rb(idx) + (xi(idx)*mu(2)*ud);
                
                %% Left & Right
                % right
                pos = this.hlp.yd <= this.hlp.yr*mu(1)/2;
                idx = this.idxmat(end,pos);
                rb(idx) = rb(idx) + (xi(idx)*mu(2)*ud);
                % left
                pos = this.hlp.yd <= this.hlp.yr*mu(1)/2;
                idx = this.idxmat(1,pos);
                rb(idx) = rb(idx) + (xi(idx)*mu(2)*ud);
                
                % Indices:
                % km3 = 14, km8 = 15, km12 = 19
                
                % x_a
                fx(1:m) = rc(2)*xi.*ya - rc(5)*xa -rc(11)*xa.*bar + rc(18)*xab + rb/this.hlp.hs;
                % y_a
                fx(m+1:2*m) = rc(1)*yi.*xan - rc(6)*ya - rc(3)*ya.*iap + rc(14)*yab;
                % x_i
                fx(2*m+1:3*m) = -rc(2)*xi.*ya - rc(9)*xi + rc(16) - rb/this.hlp.hs;
                % y_i
                fx(3*m+1:4*m) = -rc(1)*yi.*xan - rc(10)*yi + rc(17);
                % iap
                fx(4*m+1:5*m) = -(rc(3)+rc(4))*ya.*iap - rc(8)*iap + rc(15) + rc(14)*yab;
                % bar
                fx(5*m+1:6*m) = -rc(11)*xa.*bar + rc(18)*xab - rc(19)*bar + rc(19);
                % yab
                fx(6*m+1:7*m) = rc(3)*ya.*iap -(rc(14)+rc(7))*yab;
                % xab
                fx(7*m+1:8*m) = rc(11)*xa.*bar-(rc(18)+rc(13))*xab;
            else
                % Extract single functions
                xa = x(1:m,:);
                xan = bsxfun(@power,xa,mu(4,:));
                ya = x(m+1:2*m,:);
                xi = x(2*m+1:3*m,:);
                yi = x(3*m+1:4*m,:); 
                iap = x(4*m+1:5*m,:);
                bar = x(5*m+1:6*m,:);
                yab = x(6*m+1:7*m,:);
                xab = x(7*m+1:end,:);

                % Compile boundary conditions
                nd = size(x,2);
                rb = zeros(m,nd);
                
                %% Top & Bottom
                xd = repmat(this.hlp.xd,1,nd); % y distances
                % bottom
                pos = bsxfun(@lt,xd,this.hlp.xr*mu(1,:)/2);
                idx = this.idxmat(:,1);
                rb(idx,:) = pos .* (bsxfun(@times,xi(idx,:),mu(2,:).*ud));
                % top
                pos = bsxfun(@lt,xd,this.hlp.xr*mu(1,:)/2);
                idx = this.idxmat(:,end);
                rb(idx,:) = rb(idx,:) + pos .* (bsxfun(@times,xi(idx,:),mu(2,:).*ud));
                
                %% Left & Right
                yd = repmat(this.hlp.yd,1,nd);
                % right
                pos = bsxfun(@lt,yd,this.hlp.yr*mu(1,:)/2);
                idx = this.idxmat(end,:);
                rb(idx,:) = rb(idx,:) + pos .* (bsxfun(@times,xi(idx,:),mu(2,:).*ud));
                % left
                pos = bsxfun(@lt,yd,this.hlp.yr*mu(1,:)/2);
                idx = this.idxmat(1,:);
                rb(idx,:) = rb(idx,:) + pos .* (bsxfun(@times,xi(idx,:),mu(2,:).*ud));
                
                % x_a
                fx(1:m,:) = rc(2)*xi.*ya - rc(5)*xa -rc(11)*xa.*bar + rc(18)*xab + rb/this.hlp.hs;
                % y_a
                fx(m+1:2*m,:) = rc(1)*yi.*xan - rc(6)*ya - rc(3)*ya.*iap + rc(14)*yab;
                % x_i
                fx(2*m+1:3*m,:) = -rc(2)*xi.*ya - rc(9)*xi + rc(16) - rb/this.hlp.hs;
                % y_i
                fx(3*m+1:4*m,:) = -rc(1)*yi.*xan - rc(10)*yi + rc(17);
                % iap
                fx(4*m+1:5*m,:) = -(rc(3)+rc(4))*ya.*iap - rc(8)*iap + rc(15) + rc(14)*yab;
                % bar
                fx(5*m+1:6*m,:) = -rc(11)*xa.*bar + rc(18)*xab - rc(19)*bar + rc(19);
                % yab
                fx(6*m+1:7*m,:) = rc(3)*ya.*iap -(rc(14)+rc(7))*yab;
                % xab
                fx(7*m+1:8*m,:) = rc(11)*xa.*bar-(rc(18)+rc(13))*xab;
            end
            
            % If this has been projected, project back to reduced space
            if ~isempty(this.W)
                fx = this.W'*fx;
            end
        end
        
        function J = getStateJacobian(this, x, t, mu)
            
            % If this has been projected, restore full size and compute values.
            if ~isempty(this.V)
                x = this.V*x;
            end
            
            n = this.nodes;
            rc = this.System.ReacCoeff;
            
            % Boundary stuff
            bottom = this.idxmat(this.hlp.xd <= this.hlp.xr*mu(1)/2,1);
            top = this.idxmat(this.hlp.xd <= this.hlp.xr*mu(1)/2,end);
            right = this.idxmat(end,this.hlp.yd <= this.hlp.yr*mu(1)/2);
            left = this.idxmat(1,this.hlp.yd <= this.hlp.yr*mu(1)/2);
            rbxi = zeros(n,1);
            u = this.activationFun(t);
            rbxi(bottom) = mu(2)*u;
            rbxi(top) = mu(2)*u;
            rbxi(right) = mu(2)*u;
            rbxi(left) = mu(2)*u;
            rbxi = rbxi/this.hlp.hs;
            
            %% 1=x_a, 2=y_a, 3=x_i, 4=y_i, 5=iap, 6=bar, 7=yb, 8=xb
            % Add x_a dependencies
            % rc(2)*xi.*ya - rc(5)*xa -rc(11)*xa.*bar + rc(18)*xb
            i = repmat(pos(1),1,5); 
            j = pos([1:3 6 8]);
            s = [-ones(n,1)*rc(5); ... % dx_a/x_a = -k5
                rc(2)*x(pos(3));... % dx_a/y_a = k2*x_i
                rc(2)*x(pos(2)) + rbxi; %dx_a/x_i + rbxi
                -rc(11)*xa; % dx_a/bar = -rc(11)*xa
                ones(n,1)*rc(18)]; %dx_a/xb = rc(18)
            
            % Add y_a dependencies
            % rc(1)*yi.*xa - rc(6)*ya - rc(3)*ya.*iap + rc(14)*yb;
            i = [i repmat(pos(2),1,5)];
            j = [j pos([1:2 4 5 7])];
            dya_xy = mu(4)*rc(1)*x(pos(4)).*x(1:n).^(mu(4)-1);
            k1xan = rc(1)*x(1:n).^mu(4);
            s = [s; dya_xy; ... % dy_a/x_a = n*k1*y_i*x_a^(n-1)
                -ones(n,1)*rc(6);... % dy_a/y_a = -k6
                k1xan;... %dy_a/y_i = k1 * x_a^n
                -rc(3)*x(pos(2));...
                ones(n,1)*rc(14)];
            
            %% 1=x_a, 2=y_a, 3=x_i, 4=y_i, 5=iap, 6=bar, 7=yb, 8=xb
            % Add x_i dependencies
            % -rc(2)*xi.*ya - rc(9)*xi + rc(16) - rb/this.hlp.hs;
            i = [i repmat(pos(3),1,2)]; 
            j = [j pos([2 3])];
            s = [s; -rc(2)*x(pos(3));... %dx_i/y_a = -k2*x_i
                -rc(2)*x(pos(2))-rc(9)-rbxi]; %dx_i/x_i = -k2*y_a -k9 -rbxi
                
            % Add y_i dependencies HIER WEITER
            % -rc(1)*yi.*xa - rc(10)*yi + rc(17);
            i = [i repmat(pos(4),1,2)]; 
            j = [j pos([1 4])];
            s = [s; -dya_xy; ... %dy_i/x_a = -n*k1*y_i*x_a^(n-1)
                -k1xan -rc(10)]; % dy_i/y_i = -k1 * x_a^n-k10
            
            % 1=x_a, 2=y_a, 3=x_i, 4=y_i, 5=iap, 6=bar, 7=yb, 8=xb
            % Add iap dependencies
            % -(rc(3)+rc(4))*ya.*iap - rc(8)*iap + rc(15) + rc(14)*yab;
            i = [i repmat(pos(5),1,3)]; 
            j = [j pos([2 5 7])];
            % Add bar dependencies
            % -rc(11)*xa.*bar + rc(18)*xab - rc(19)*bar + rc(19);
            i = [i repmat(pos(6),1,3)]; 
            j = [j pos([1 6 8])];
            % Add bar dependencies
            % rc(3)*ya.*iap -(rc(14)+rc(7))*yab;
            i = [i repmat(pos(7),1,3)]; 
            j = [j pos([2 5 7])];
            % Add bar dependencies
            % rc(11)*xa.*bar-(rc(18)+rc(13))*xab;
            i = [i repmat(pos(8),1,3)]; 
            j = [j pos([1 6 8])];
            
            n = this.fDim;
            if ~isempty(this.V)
                n = size(this.V,1);
            end
            m = this.xDim;
            if ~isempty(this.W)
                n = size(this.W,1);
            end
            J = sparse(i,j,s,n,m);
            
            function idx = pos(nr)
                idx = [];
                for k=1:length(nr)
                    idx = [idx (nr(k)-1)*n+1:nr(k)*n];%#ok
                end
            end
        end
    end
    
    methods(Access=protected)
        function fxj = evaluateComponents(this, pts, ends, ~, ~, X, t, mu)
            % The vector embedding results from the fixed ordering of the full 4*m-vector into
            % the components x_a, y_a, x_i, y_i
            %
            % Parameters:
            % pts: The components of `\vf` for which derivatives are required @type rowvec<integer>
            % ends: At the `i`-th entry it contains the last position in the `\vx` vector that
            % indicates an input value relevant for the `i`-th point evaluation, i.e.
            % `f_i(\vx) = f_i(\vx(ends(i-1){:}ends(i)));` @type rowvec<integer>
            % X: A matrix `\vX` with the state space locations `\vx_i` in its columns @type
            % matrix<double>
            % mu: The corresponding parameters `\mu_i` for each state `\vx_i`, as column matrix
            % @type matrix<double>
            %
            % Return values:
            % fxj: A matrix with pts-many component function evaluations `f_i(\vx)` as rows and as
            % many columns as `\vX` had.
            m = this.nodes;
            nd = size(X,2);
            fxj = zeros(length(pts),nd);
            rc = this.System.ReacCoeff;
            if nd > 1
                bottom = bsxfun(@lt,this.hlp.xd,this.hlp.xr*mu(1,:)/2);
                top = bsxfun(@lt,this.hlp.xd,this.hlp.xr*mu(1,:)/2);
                right = bsxfun(@lt,this.hlp.yd,this.hlp.yr*mu(1,:)/2);
                left = bsxfun(@lt,this.hlp.yd,this.hlp.yr*mu(1,:)/2);
            else
                bottom = this.hlp.xd <= this.hlp.xr*mu(1,:)/2;
                top = this.hlp.xd <= this.hlp.xr*mu(1,:)/2;
                right = this.hlp.yd <= this.hlp.yr*mu(1,:)/2;
                left = this.hlp.yd <= this.hlp.yr*mu(1,:)/2;
            end
            % Get matrix indices
            J2 = mod(pts,m);
            J2(J2==0) = m;
            
            %[row2, col2] = ind2sub([this.hlp.d1 this.hlp.d2],J2);
            % ind2sub direct replacement!
            row = rem(J2-1, this.hlp.d1)+1;
            col = (J2-row)/this.hlp.d1+1;
            u = this.activationFun(t, mu);
            for idx=1:length(pts)
                j = pts(idx);
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
                    % rc(1)*xi.*ya - rc(3)*xa + rb;
                    fj = rc(1)*x(3,:).*x(2,:) - rc(3)*x(1,:);
                    
                    % Boundary conditions
                    % Bottom
                    if col(idx) == 1
                        fj = fj + bottom(row(idx),:) .* (x(3,:).*mu(2,:).*u/this.hlp.hs);
                    end
                    % Top
                    if col(idx) == this.hlp.d2
                        fj = fj + top(row(idx),:) .* (x(3,:).*mu(2,:).*u/this.hlp.hs);
                    end
                    % Right
                    if row(idx) == this.hlp.d1
                        fj = fj + right(col(idx),:) .* (x(3,:).*mu(2,:).*u/this.hlp.hs);
                    end
                    % Left
                    if row(idx) == 1
                        fj = fj + left(col(idx),:) .* (x(3,:).*mu(2,:).*u/this.hlp.hs);
                    end
                    
                % Y_a
                elseif m < j && j <= 2*m
                    % 1=x_a, 2=y_a {, x_i}, 3=y_i
                    % rc(2)*yi.*xa^mu4 - rc(4)*ya;
                    fj = rc(2)*x(3,:).*x(1,:).^mu(4,:) - rc(4)*x(2,:);
                    
                % X_i
                elseif 2*m < j && j <= 3*m
                    % {x_a,} 1=y_a, 2=x_i {, y_i}
                    % -rc(1)*xi.*ya - rc(5)*xi + rc(7) - rb;
                    fj = -rc(1)*x(2,:).*x(1,:) - rc(5)*x(2,:) + rc(7);
                    
                    % Boundary conditions
                    % Bottom
                    if col(idx) == 1
                        fj = fj - bottom(row(idx),:) .* (x(2,:).*mu(2,:).*u/this.hlp.hs);
                    end
                    % Top
                    if col(idx) == this.hlp.d2
                        fj = fj - top(row(idx),:) .* (x(2,:).*mu(2,:).*u/this.hlp.hs);
                    end
                    % Right
                    if row(idx) == this.hlp.d1
                        fj = fj - right(col(idx),:) .* (x(2,:).*mu(2,:).*u/this.hlp.hs);
                    end
                    % Left
                    if row(idx) == 1
                        fj = fj - left(col(idx),:) .* (x(2,:).*mu(2,:).*u/this.hlp.hs);
                    end
                    
                % Y_i
                else
                    % 1=x_a, {y_a, x_i}, 2=y_i
                    % -rc(2)*yi.*xa^mu4 - rc(6)*yi + rc(8);
                    fj = -rc(2)*x(2,:).*x(1,:).^mu(4,:) - rc(6)*x(2,:) + rc(8);
                end
                fxj(idx,:) = fj;
            end
        end
        
        function dfx = evaluateComponentPartialDerivatives(this, pts, ends, ~, deriv, ~, X, ~, mu, ~)
            % See dscomponents.ACompEvalCoreFun for more details.
            %
            % Parameters:
            % pts: The components of `f` for which derivatives are required @type
            % rowvec<integer>
            % ends: At the `i`-th entry it contains the last position in the `\vx` vector that
            % indicates an input value relevant for the `i`-th point evaluation, i.e.
            % `f_i(\vx) = f_i(\vx(ends(i-1){:}ends(i)));` @type rowvec<integer>
            % deriv: The indices within `\vx` that derivatives are required for.
            % @type rowvec<integer>
            % X: The state space location `\vx` @type colvec<double>
            % mu: The corresponding parameter `\mu` for the state `\vx` @type colvec<double>
            %
            % Return values:
            % dfx: A column vector with 'numel(deriv)' rows containing the derivatives at all
            % specified pts i with respect to the coordinates given by 'idx(ends(i-1):ends(i))'
            error('Activation fun not yet implemented correctly');
            m = this.nodes;
            rc = this.System.ReacCoeff;
            nd = size(X,2);
            
            %% Boundary stuff
            % Get matrix indices of points
            pts2 = mod(pts,m);
            pts2(pts2==0) = m;
            row = rem(pts2-1, this.hlp.d1)+1;
            col = (pts2-row)/this.hlp.d1+1;
            if nd > 1
                bottom = bsxfun(@lt,this.hlp.xd,this.hlp.xr*mu(1)/2);
                top = bsxfun(@lt,this.hlp.xd,this.hlp.xr*mu(1)/2);
                right = bsxfun(@lt,this.hlp.yd,this.hlp.yr*mu(1)/2);
                left = bsxfun(@lt,this.hlp.yd,this.hlp.yr*mu(1)/2);
            else
                bottom = this.hlp.xd <= this.hlp.xr*mu(1)/2;
                top = this.hlp.xd <= this.hlp.xr*mu(1)/2;
                right = this.hlp.yd <= this.hlp.yr*mu(1)/2;
                left = this.hlp.yd <= this.hlp.yr*mu(1)/2;
            end
            
            %% Derivative info per point
            der = false(size(X,1),1);
            der(deriv) = true;
           
            %% Main loop
            dfx = zeros(size(deriv,1),nd);
            curpos = 1;
            for idx=1:length(pts)
                i = pts(idx);
                if idx == 1
                    st = 0;
                else
                    st = ends(idx-1);
                end
                % Select the elements of x that are effectively used in f
                elem = (st+1):ends(idx);
                x = X(elem,:);
                d = der(elem);
                
                % X_a
                if i <= m
                    % 1=x_a, 2=y_a, 3=x_i {, y_i}
                    if d(1) % dx_a/x_a = -mu3
                        dfx(curpos,:) = -rc(3);
                        curpos = curpos + 1;
                    end
                    if d(2) % dx_a/y_a = mu1*x_i
                        dfx(curpos,:) = rc(1)*x(3,:);
                        curpos = curpos + 1;
                    end
                    if d(3) % dx_a/x_i = mu1*y_a + rb'
                        dfx(curpos,:) = rc(1)*x(2,:);
                        if col(idx) == 1 % Bottom
                            dfx(curpos,:) = dfx(curpos,:) + bottom(row(idx),:) .* mu(2,:)/this.hlp.hs;
                        end
                        if col(idx) == this.hlp.d2 % Top
                            dfx(curpos,:) = dfx(curpos,:) + top(row(idx),:) .* mu(2,:)/this.hlp.hs;
                        end
                        if row(idx) == this.hlp.d1 % Right
                            dfx(curpos,:) = dfx(curpos,:) + right(col(idx),:) .* mu(2,:)/this.hlp.hs;
                        end
                        if row(idx) == 1 % Left
                            dfx(curpos,:) = dfx(curpos,:) + left(col(idx),:) .* mu(2,:)/this.hlp.hs;
                        end
                        curpos = curpos + 1;
                    end
                    
                % Y_a
                elseif m < i && i <= 2*m
                    % 1=x_a, 2=y_a {, x_i}, 3=y_i
                    if d(1) % dy_a/x_a = n*mu2*y_i*x_a^(n-1)
                        dfx(curpos,:) = mu(4,:)*rc(2)*x(3,:).*x(1,:).^(mu(4,:)-1);
                        curpos = curpos + 1;
                    end
                    if d(2) % dy_a/y_a = -mu4
                        dfx(curpos,:) = -rc(4);
                        curpos = curpos + 1;
                    end
                    if d(3) % dx_a/x_a = -mu3
                        dfx(curpos,:) = rc(2)*x(1,:).^mu(4,:);
                        curpos = curpos + 1;
                    end
                    
                % X_i
                elseif 2*m < i && i <= 3*m
                    % {x_a,} 1=y_a, 2=x_i {, y_i}
                    if d(1) %dx_i/y_a = -mu1*x_i
                        dfx(curpos,:) = -rc(1)*x(2,:);
                        curpos = curpos + 1;
                    end
                    if d(2) %dx_i/x_i = -mu1*y_a -mu5 -rbxi
                        dfx(curpos,:) = -rc(1)*x(1,:)-rc(5);
                        % Boundary conditions
                        if col(idx) == 1 % Bottom
                            dfx(curpos,:) = dfx(curpos,:) - bottom(row(idx),:) .* mu(2,:)/this.hlp.hs;
                        end
                        if col(idx) == this.hlp.d2 % Top
                            dfx(curpos,:) = dfx(curpos,:) - top(row(idx),:) .* mu(2,:)/this.hlp.hs;
                        end
                        if row(idx) == this.hlp.d1 % Right
                            dfx(curpos,:) = dfx(curpos,:) - right(col(idx),:) .* mu(2,:)/this.hlp.hs;
                        end
                        if row(idx) == 1 % Left
                            dfx(curpos,:) = dfx(curpos,:) - left(col(idx),:) .* mu(2,:)/this.hlp.hs;
                        end
                        curpos = curpos + 1;
                    end
                    
                % Y_i
                else
                    % 1=x_a, {y_a, x_i}, 2=y_i
                    if d(1) %dy_i/x_a = -n*mu2*y_i*x_a^(n-1)
                        dfx(curpos,:) = -mu(4,:)*rc(2)*x(2,:).*x(1,:).^(mu(4,:)-1);
                        curpos = curpos + 1;
                    end
                    if d(2) % dy_i/y_i = -mu2 * x_a^n-mu6
                        dfx(curpos,:) = -rc(2)*x(1,:).^mu(4,:) - rc(6);
                        curpos = curpos + 1;
                    end
                end
            end
        end
    end
end

