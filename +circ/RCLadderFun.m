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
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    methods
        function this = RCLadderFun(sys, dim)
            this = this@dscomponents.ACoreFun(sys);
            this.CustomProjection = true;
            this.TimeDependent = false;
            this.xDim = dim;
            this.fDim = dim;
        end
        
        function fx = evaluateCoreFun(~, v, ~)
            vd = exp(40*(v(1:end-1,:)-v(2:end,:)));
            vdf = [vd; ones(1,size(v,2))];
            vdb = [-exp(40*v(1,:)) + 2; vd];
            
            fx = vdb - vdf;
        end
        
        function fx = evaluateMulti(~, v, ~, ~)
            vd = exp(40*(v(1:end-1,:)-v(2:end,:)));
            vdf = [vd; ones(1,size(v,2))];
            vdb = [-exp(40*v(1,:)) + 2; vd];
            
            fx = vdb - vdf;
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