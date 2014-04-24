classdef InhibitCoreFun2D < models.pcdi.BaseCoreFun
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
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    methods
        
        function this = InhibitCoreFun2D(dynsys)
            this = this@models.pcdi.BaseCoreFun(dynsys);
        end
           
        function copy = clone(this)
            % System already copied in constructor (see below)
            copy = models.pcdi.InhibitCoreFun2D(this.System);
            
            % Call superclass method
            copy = clone@models.pcdi.BaseCoreFun(this, copy);
        end
        
        function newSysDimension(this)
            % Create diffusion matrix
            s = this.System;
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
            this.hlp.xd = abs(((1:this.hlp.d1)-1)*this.System.h-.5*this.hlp.xr)'; % x distances from middle
            this.hlp.yd = abs(((1:this.hlp.d2)-1)*this.System.h-.5*this.hlp.yr)'; % y distances from middle
            
            % Build sparsity pattern
            n = this.nodes;
            
            %% 1=x_a, 2=y_a, 3=x_i, 4=y_i, 5=iap, 6=bar, 7=yb, 8=xb
            % Add x_a dependencies
            % rc(2)*xi.*ya - rc(5)*xa -rc(11)*xa.*bar + rc(18)*xb
            i = repmat(this.nodepos(1),1,5); 
            j = this.nodepos([1:3 6 8]);
            % Add y_a dependencies
            % rc(1)*yi.*xa - rc(6)*ya - rc(3)*ya.*iap + rc(14)*yb;
            i = [i repmat(this.nodepos(2),1,5)];
            j = [j this.nodepos([1:2 4 5 7])];
            % Add x_i dependencies
            % -rc(2)*xi.*ya - rc(9)*xi + rc(16) - rb/this.hlp.hs;
            i = [i repmat(this.nodepos(3),1,2)]; 
            j = [j this.nodepos([2 3])];
            % Add y_i dependencies
            % -rc(1)*yi.*xa - rc(10)*yi + rc(17);
            i = [i repmat(this.nodepos(4),1,2)]; 
            j = [j this.nodepos([1 4])];
            
            %% 1=x_a, 2=y_a, 3=x_i, 4=y_i, 5=iap, 6=bar, 7=yb, 8=xb
            % Add iap dependencies
            % -(rc(3)+rc(4))*ya.*iap - rc(8)*iap + rc(15) + rc(14)*yab;
            i = [i repmat(this.nodepos(5),1,3)]; 
            j = [j this.nodepos([2 5 7])];
            % Add bar dependencies
            % -rc(11)*xa.*bar + rc(18)*xab - rc(19)*bar + rc(19);
            i = [i repmat(this.nodepos(6),1,3)]; 
            j = [j this.nodepos([1 6 8])];
            % Add bar dependencies
            % rc(3)*ya.*iap -(rc(14)+rc(7))*yab;
            i = [i repmat(this.nodepos(7),1,3)]; 
            j = [j this.nodepos([2 5 7])];
            % Add bar dependencies
            % rc(11)*xa.*bar-(rc(18)+rc(13))*xab;
            i = [i repmat(this.nodepos(8),1,3)]; 
            j = [j this.nodepos([1 6 8])];
            this.Ji = i;
            this.Jj = j;
            this.xDim = 8*n;
            this.fDim = 8*n;
            this.JSparsityPattern = sparse(i,j,ones(length(i),1),this.fDim,this.xDim);
        end
        
        function fx = evaluate(this, x, t)
            % If this has been projected, restore full size and compute values.
            if ~isempty(this.V)
                x = this.V*x;
            end
            
            % Allocate result vector
            fx = zeros(size(x));

            m = this.nodes;
            mu = this.mu;
            diff = this.System.CurCXMU;
            
            % Compute activation fun values (it's set up in scaled times!)
            ud = this.activationFun(t,mu);
            rc = this.System.ReacCoeff;
            
            % Extract single functions
            xa = x(1:m);
            xan = xa.^mu(4);
            ya = x(m+1:2*m);
            xi = x(2*m+1:3*m);
            yi = x(3*m+1:4*m);
            iap = x(4*m+1:5*m);
            bar = x(5*m+1:6*m);
            yb = x(6*m+1:7*m);
            xb = x(7*m+1:end);

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
            fx(1:m) = rc(2)*xi.*ya - rc(5)*xa -rc(11)*xa.*bar + rc(18)*xb + rb/this.hlp.hs;
            % y_a
            fx(m+1:2*m) = rc(1)*yi.*xan - rc(6)*ya - rc(3)*ya.*iap + rc(14)*yb;
            % x_i
            fx(2*m+1:3*m) = -rc(2)*xi.*ya - rc(9)*xi + rc(16)*diff - rb/this.hlp.hs;
            % y_i
            fx(3*m+1:4*m) = -rc(1)*yi.*xan - rc(10)*yi + rc(17)*diff;
            % iap
            fx(4*m+1:5*m) = -(rc(3)+rc(4))*ya.*iap - rc(8)*iap + rc(15)*diff + rc(14)*yb;
            % bar
            fx(5*m+1:6*m) = -rc(11)*xa.*bar + rc(18)*xb - rc(19)*bar + rc(19)*diff;
            % yab
            fx(6*m+1:7*m) = rc(3)*ya.*iap -(rc(14)+rc(7))*yb;
            % xab
            fx(7*m+1:8*m) = rc(11)*xa.*bar-(rc(18)+rc(13))*xb;
            
            % If this has been projected, project back to reduced space
            if ~isempty(this.W)
                fx = this.W'*fx;
            end
        end
        
        function fx = evaluateMulti(this, x, t, mu)
            error('not usable due to spatially dependent production rates.');
            
            % Allocate result vector
            fx = zeros(size(x));

            m = this.nodes;
            
            % Compute activation fun values (it's set up in scaled times!)
            ud = this.activationFun(t,mu);
            rc = this.System.ReacCoeff;
            
            % Extract single functions
            xa = x(1:m,:);
            xan = bsxfun(@power,xa,mu(4,:));
            ya = x(m+1:2*m,:);
            xi = x(2*m+1:3*m,:);
            yi = x(3*m+1:4*m,:);
            iap = x(4*m+1:5*m,:);
            bar = x(5*m+1:6*m,:);
            yb = x(6*m+1:7*m,:);
            xb = x(7*m+1:end,:);
            
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
            fx(1:m,:) = rc(2)*xi.*ya - rc(5)*xa -rc(11)*xa.*bar + rc(18)*xb + rb/this.hlp.hs;
            % y_a
            fx(m+1:2*m,:) = rc(1)*yi.*xan - rc(6)*ya - rc(3)*ya.*iap + rc(14)*yb;
            % x_i
            fx(2*m+1:3*m,:) = -rc(2)*xi.*ya - rc(9)*xi + rc(16) - rb/this.hlp.hs;
            % y_i
            fx(3*m+1:4*m,:) = -rc(1)*yi.*xan - rc(10)*yi + rc(17);
            % iap
            fx(4*m+1:5*m,:) = -(rc(3)+rc(4))*ya.*iap - rc(8)*iap + rc(15) + rc(14)*yb;
            % bar
            fx(5*m+1:6*m,:) = -rc(11)*xa.*bar + rc(18)*xb - rc(19)*bar + rc(19);
            % yab
            fx(6*m+1:7*m,:) = rc(3)*ya.*iap -(rc(14)+rc(7))*yb;
            % xab
            fx(7*m+1:8*m,:) = rc(11)*xa.*bar-(rc(18)+rc(13))*xb;
        end
        
        function J = getStateJacobian(this, x, t)
            
            % If this has been projected, restore full size and compute values.
            if ~isempty(this.V)
                x = this.V*x;
            end
            
            n = this.nodes;
            mu = this.mu;
            rc = this.System.ReacCoeff;
            
            % Boundary stuff
            bottom = this.idxmat(this.hlp.xd <= this.hlp.xr*mu(1)/2,1);
            top = this.idxmat(this.hlp.xd <= this.hlp.xr*mu(1)/2,end);
            right = this.idxmat(end,this.hlp.yd <= this.hlp.yr*mu(1)/2);
            left = this.idxmat(1,this.hlp.yd <= this.hlp.yr*mu(1)/2);
            rbxi = zeros(n,1);
            u = this.activationFun(t,mu);
            rbxi(bottom) = mu(2)*u;
            rbxi(top) = mu(2)*u;
            rbxi(right) = mu(2)*u;
            rbxi(left) = mu(2)*u;
            rbxi = rbxi/this.hlp.hs;
            
            % yb, xb values not needed in jacobian!
            hlp = reshape(1:6*n',[],6)';
            o = ones(n,1);
            xa = x(hlp(1,:)); ya = x(hlp(2,:));
            xi = x(hlp(3,:)); yi = x(hlp(4,:));
            iap = x(hlp(5,:)); bar = x(hlp(6,:));
            
            %% 1=x_a, 2=y_a, 3=x_i, 4=y_i, 5=iap, 6=bar, 7=yb, 8=xb
            % dxa = rc(2)*xi.*ya - rc(5)*xa -rc(11)*xa.*bar + rc(18)*xb
            s = [-o*rc(5)-rc(11)*bar; ... % dx_a/x_a
                rc(2)*xi;... % dxa/ya
                rc(2)*ya + rbxi; %dxa/xi
                -rc(11)*xa; % dxa/bar
                o*rc(18)]; %dxa/xb
            % dya = rc(1)*yi.*xan - rc(6)*ya - rc(3)*ya.*iap + rc(14)*yb;
            nyixan_1 = mu(4)*rc(1)*yi.*xa.^(mu(4)-1);
            k1xan = rc(1)*xa.^mu(4);
            s = [s; nyixan_1; ... % dya/xa
                -o*rc(6) - rc(3)*iap;... % dy_a/y_a
                k1xan;... %dy_a/y_i
                -rc(3)*xa;... %dya/
                o*rc(14)];
            %% 1=x_a, 2=y_a, 3=x_i, 4=y_i, 5=iap, 6=bar, 7=yb, 8=xb
            % Add x_i dependencies
            % dxi = -rc(2)*xi.*ya - rc(9)*xi + rc(16) - rb/this.hlp.hs;
            s = [s; -rc(2)*xi;... %dxi/ya
                -rc(2)*xa-rc(9)-rbxi]; %dxi/xi
            % -rc(1)*yi.*xan - rc(10)*yi + rc(17);
            s = [s; -nyixan_1; ... %dyi/xa
                -k1xan-o*rc(10)]; % dyi/yi
            
            %% 1=x_a, 2=y_a, 3=x_i, 4=y_i, 5=iap, 6=bar, 7=yb, 8=xb
            % Add iap dependencies
            % diap = -(rc(3)+rc(4))*ya.*iap - rc(8)*iap + rc(15) + rc(14)*yb;
            s = [s; -(rc(3)+rc(4))*iap; ... %diap/ya
                -(rc(3)+rc(4))*ya-rc(8); % diap/iap
                o*rc(14)]; % diap/yb
            % dbar = -rc(11)*xa.*bar + rc(18)*xb - rc(19)*bar + rc(19);
            s = [s; -rc(11)*bar; ... %dbar/xa
                -rc(11)*xa-rc(8); % dbar/bar
                o*rc(18)]; % diap/xb
            %% 1=x_a, 2=y_a, 3=x_i, 4=y_i, 5=iap, 6=bar, 7=yb, 8=xb
            % dyb = rc(3)*ya.*iap -(rc(14)+rc(7))*yb;
            s = [s; rc(3)*iap; ... %dyb/ya
                rc(3)*ya; % dyb/iap
                -o*(rc(14)+rc(7))]; % dyb/yb
            % dxb = rc(11)*xa.*bar-(rc(18)+rc(13))*xb;
            s = [s; rc(11)*bar; ... %dxb/xa
                rc(11)*xa; % dxb/bar
                -o*(rc(18)+rc(13))]; % dxp/xb
            
            n = this.fDim;
            if ~isempty(this.V)
                n = size(this.V,1);
            end
            m = this.xDim;
            if ~isempty(this.W)
                n = size(this.W,1);
            end
            J = sparse(this.Ji,this.Jj,s,n,m);
        end
    end
    
    methods(Access=protected)
        function fxj = evaluateComponents(this, pts, ends, globidx, ~, X, t, mu)
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
            diff = this.System.CurCXMU;
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
                
                %% X_a
                % General: 1=x_a, 2=y_a, 3=x_i, 4=y_i, 5=iap, 6=bar, 7=yb, 8=xb
                % Sel: 1=xa, 2=ya, 3=xi, 4=bar, 5=xb
                % rc(2)*xi.*ya - rc(5)*xa -rc(11)*xa.*bar + rc(18)*xb + rb/this.hlp.hs;
                if j <= m
                    fj = rc(2)*x(3,:).*x(2,:) - rc(5)*x(1,:) -rc(11)*x(1,:).*x(4,:) + rc(18)*x(5,:);
                    
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
                    
                %% Y_a
                % General: 1=x_a, 2=y_a, 3=x_i, 4=y_i, 5=iap, 6=bar, 7=yb, 8=xb
                % Sel: 1=x_a, 2=y_a, 3=y_i, 4=iap, 5=yb
                % rc(1)*yi.*xan - rc(6)*ya - rc(3)*ya.*iap + rc(14)*yb;
                elseif m < j && j <= 2*m
                    fj = rc(1)*x(3,:).*x(1,:).^mu(4,:) - rc(6)*x(2,:) - rc(3)*x(2,:).*x(4,:) + rc(14)*x(5,:);
                    
                %% X_i
                % General: 1=x_a, 2=y_a, 3=x_i, 4=y_i, 5=iap, 6=bar, 7=yb, 8=xb
                % Sel: % 1=ya, 2=xi
                % -rc(2)*xi.*ya - rc(9)*xi + rc(16) - rb/this.hlp.hs;
                elseif 2*m < j && j <= 3*m
                    fj = -rc(2)*x(2,:).*x(1,:) - rc(9)*x(2,:) + rc(16);
                    
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
                    
                %% Y_i
                % General: 1=x_a, 2=y_a, 3=x_i, 4=y_i, 5=iap, 6=bar, 7=yb, 8=xb
                % Sel: 1=x_a, 2=y_i
                % -rc(1)*yi.*xan - rc(10)*yi + rc(17);
                elseif 3*m < j && j <= 4*m
                    fj = -rc(1)*x(2,:).*x(1,:).^mu(4,:) - rc(10)*x(2,:) + rc(17);
                    
                %% IAP
                % General: 1=x_a, 2=y_a, 3=x_i, 4=y_i, 5=iap, 6=bar, 7=yb, 8=xb
                % Sel: 1=ya, 2=iap, 3=yb
                % -(rc(3)+rc(4))*ya.*iap - rc(8)*iap + rc(15) + rc(14)*yb;
                elseif 4*m < j && j <= 5*m
                    fj = -(rc(3)+rc(4))*x(1,:).*x(2,:) - rc(8)*x(2,:) + rc(15) + rc(14)*x(3,:);
                    
                %% BAR
                % General: 1=x_a, 2=y_a, 3=x_i, 4=y_i, 5=iap, 6=bar, 7=yb, 8=xb
                % Sel: 1=xa, 2=bar, 3=xb
                % -rc(11)*xa.*bar + rc(18)*xb - rc(19)*bar + rc(19);
                elseif 5*m < j && j <= 6*m
                    fj = -rc(11)*x(1,:).*x(2,:) + rc(18)*x(3,:) - rc(19)*x(2,:) + rc(19);
                    
                %% YB
                % General: 1=x_a, 2=y_a, 3=x_i, 4=y_i, 5=iap, 6=bar, 7=yb, 8=xb
                % Sel: 1=ya, 2=iap, 3=yb
                % rc(3)*ya.*iap -(rc(14)+rc(7))*yb;
                elseif 6*m < j && j <= 7*m
                    fj = rc(3)*x(1,:).*x(2,:) -(rc(14)+rc(7))*x(3,:);
                
                %% XB
                % General: 1=x_a, 2=y_a, 3=x_i, 4=y_i, 5=iap, 6=bar, 7=yb, 8=xb
                % Sel: 1=xa, 2=bar, 3=xb
                % rc(11)*xa.*bar-(rc(18)+rc(13))*xb;
                else
                    fj = rc(11)*x(1,:).*x(2,:)-(rc(18)+rc(13))*x(3,:);
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

