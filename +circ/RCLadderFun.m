classdef RCLadderFun < dscomponents.ACoreFun
% RCLadderFun: 
%
%
%
% @author Daniel Wirtz @date 2011-04-29
%
% @new{0,3,dw,2011-04-29} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods
        function this = RCLadderFun
            this = this@dscomponents.ACoreFun;
            
            this.MultiArgumentEvaluations = true;
            this.CustomProjection = true;
        end
        
        function fx = evaluateCoreFun(this, v, t, mu)%#ok
            
            d = size(v,1);
            e = ones(d,1);
            A = spdiags([e -2*e e],-1:1,d,d);
            
            sv = sign(v);
            sv(sv == 0) = 1;
            fx = A*v - sv.*v.^2;
            
%             vd = v(1:end-1,:) - v(2:end,:);
%             %vd = v(2:end,:)-v(1:end-1,:);
%             
%             if isempty(this.V)
%                 fx = -(exp(40*vd)-1+vd);
%                 % First one gets added something
%                 %fx(1,:) = fx(1,:) + exp(40*v(1,:))-1+v(1,:);
%             else
%                 fx = -this.W'*(exp(40*this.V*vd) +1 -this.W'*(this.V*vd));
%             end
%             % Last ones is negative of one before
%             fx(end+1,:) = -fx(end,:);
        end
        
        function target = project(this, V, W)
            target = this.clone;
            target = project@dscomponents.ACoreFun(this, V, W, target); 
        end
        
        function copy = clone(this)
            copy = models.RCLadderFun(this.Model);
            copy = clone@dscomponents.ACoreFun(this, copy);
        end
    end
    
end