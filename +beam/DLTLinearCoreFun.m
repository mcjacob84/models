classdef DLTLinearCoreFun < models.beam.DLTBaseCoreFun & dscomponents.LinearCoreFun
% DLTLinearCoreFun: 
%
%
%
% @author Daniel Wirtz @date 2011-12-07
%
% @new{0,6,dw,2011-12-07} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods
        function this = DLTLinearCoreFun(system)
            this = this@models.beam.DLTBaseCoreFun(system);
        end
        
        function prepareConstants(this, ~)
            prepareConstants@models.beam.DLTBaseCoreFun(this);
            % Assemble matrix (K is constant for this case)
            this.A = -this.B_big;
            % Assign offset
            this.b = this.f_big;
        end
    end
    
    methods(Access=protected)
        function f_eff = getf_eff(~, f, K, u)
            f_eff = f - K*u;% + M_T*u;
        end
        
        function B = getB_big(~, K, C)
            B = [zeros(size(K)), -eye(size(K)); K, C];
        end
    end
    
end