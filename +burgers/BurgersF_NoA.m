classdef BurgersF_NoA < dscomponents.ACompEvalCoreFun
% BurgersF: 
%
%
%
% @author Daniel Wirtz @date 2012-04-24
%
% @new{0,6,dw,2012-04-24} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        Ax;
        System;
    end
    
    methods
        function this = BurgersF_NoA(sys)
            this = this@dscomponents.ACompEvalCoreFun;
            this.System = sys;
            this.MultiArgumentEvaluations = true;
            this.CustomProjection = false;
            this.TimeDependent = false;
        end
        
        function fx = evaluateCoreFun(this, x, ~, ~)
            fx = - x.*(this.Ax*x);
        end
                
        function J = getStateJacobian(this, x, ~, ~)
            hlp = bsxfun(@times,this.Ax,x);
            J = -hlp - spdiags(this.Ax*x,0,size(x,1),size(x,1));
        end
        
        function newDim(this)
            m = this.System.Model;
            n = m.Dimension;
            dx = (m.Omega(2) - m.Omega(1))/(n+1);
            e = ones(n,1);
            d1 = e/(2*dx);
            this.Ax = spdiags([-d1 0*d1  d1], -1:1, n, n);
            this.JSparsityPattern = spdiags([e e  e], -1:1, n, n);
            this.xDim = n;
            this.fDim = n;
        end
        
%         function target = project(this, V, W)
%             target = this.clone;
%             target = project@dscomponents.ACoreFun()
%         end
        
        function copy = clone(this)
            copy = clone@dscomponents.ACompEvalCoreFun(this, models.burgers.BurgersF_NoA(this.System));
            copy.Ax = this.Ax;
        end
    end
    
    methods(Access=protected)
        function fxj = evaluateComponents(this, pts, ends, argidx, self, X, ~, ~)
            % Evaluates the burgers nonlinearity pointwise.
            %
            % Parameters:
            % pts: The components of `\vf` for which derivatives are required @type rowvec<integer>
            % ends: At the `i`-th entry it contains the last position in the `\vx` vector that
            % indicates an input value relevant for the `i`-th point evaluation, i.e.
            % `f_i(\vx) = f_i(\vx(ends(i-1){:}ends(i)));` @type rowvec<integer>
            % argidx: The indices of `\vx`-entries in the global `\vx` vector w.r.t the `i`-th
            % point, e.g. `xglobal(i-1:i+1) = \vx(ends(i-1):ends(i))` @type rowvec<integer>
            % self: The positions in the `\vx` vector that correspond to the `i`-th output
            % dimension, if applicable (usually `f_i` depends on `x_i`, but not necessarily)
            % @type rowvec<integer>
            % X: A matrix `\vX` with the state space locations `\vx_i` in its columns @type
            % matrix<double>
            %
            % Return values:
            % fxj: A matrix with pts-many component function evaluations `f_i(\vx)` as rows and as
            % many columns as `\vX` had.
            fxj = zeros(length(pts),size(X,2));
            for idx=1:length(pts)
                pt = pts(idx);
                if idx == 1
                    st = 0;
                else
                    st = ends(idx-1);
                end
                % Select the elements of x that are effectively used in f
                xidx = (st+1):ends(idx);
                x = X(xidx,:);
                fxj(idx,:) = - x(self(xidx),:) .* (this.Ax(pt,argidx(xidx))*x);
            end
        end
    end
end