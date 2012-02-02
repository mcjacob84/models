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
            this.TimeDependent = false;
        end
        
        function fx = evaluateCoreFun(~, v, ~, ~)
            d = size(v,1);
            e = ones(d,1);
            A = spdiags([e -2*e e],-1:1,d,d);
            A(1,:) = -A(1,:);
            A(end,end) = -1;
            
            vd = exp(40*(v(1:end-1,:)-v(2:end,:)));
            vdf = [vd; ones(1,size(v,2))];
            vdb = [-exp(40*v(1,:)) + 2; vd];
            
            fx = vdb - vdf + A*v;
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