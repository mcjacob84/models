classdef BurgersF < dscomponents.ACompEvalCoreFun
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
        A;
        Ax;
        System;
    end
    
    methods
        function this = BurgersF(sys)
            this = this@dscomponents.ACompEvalCoreFun;
            this.System = sys;
            this.MultiArgumentEvaluations = true;
            this.CustomProjection = false;
            this.TimeDependent = false;
        end
        
        function fx = evaluateCoreFun(this, x, ~, mu)
            fx = bsxfun(@times,this.A*x,mu(1,:)) - x.*(this.Ax*x);
        end
                
        function fxj = evaluateComponents(this, pts, ends, argidx, self, X, ~, mu)
            % Evaluates the burgers nonlinearity pointwise.
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
                fxj(idx,:) = mu(1,:).*(this.A(pt,argidx(xidx))*x) - x(self(xidx),:).*(this.Ax(pt,argidx(xidx))*x);
            end
        end
        
        function J = getStateJacobian(this, x, ~, mu)
            hlp = bsxfun(@times,this.Ax,x);
            hlp = hlp + spdiags(this.Ax*x,0,size(x,1),size(x,1));
            J = mu(1)*this.A - hlp;
        end
        
        function newDim(this)
            m = this.System.Model;
            n = m.Dimension;
            dx = (m.Omega(2) - m.Omega(1))/(n+1);
            e = ones(n,1);
            d1 = e/(2*dx);
            d2 = e/(dx^2);
            this.A = spdiags([d2 -2*d2  d2], -1:1, n, n); 
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
            copy = clone@dscomponents.ACompEvalCoreFun(this, models.burgers.BurgersF(this.System));
            copy.A = this.A;
            copy.Ax = this.Ax;
        end
    end
end