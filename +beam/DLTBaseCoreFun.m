classdef DLTBaseCoreFun < dscomponents.ACoreFun
% DynLinTimoshenkoCoreFun: 
%
%
%
% @author Daniel Wirtz @date 2011-09-20
%
% @new{0,5,dw,2011-09-20} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetAccess=private, GetAccess=protected)
        sys;
    end
    
    properties(Access=protected)
        f_const_big;
    end
    
    methods
        function this = DLTBaseCoreFun(sys)
            this = this@dscomponents.ACoreFun;
            this.TimeDependent = false;
            this.sys = sys;
            this.initialize;
        end

        function initialize(this)
            % do constant preps, i.e. assign A
            % Massenmatrix (u und T)
            m = this.sys.Model;
            
            %% Jacobian matrix sparsity pattern
            % Create a fake big system matrix here to obtain the sparsity
            % pattern.
            s = length(m.free);
            K = this.sys.K0(m.free,m.free);
            null = sparse(s,s);
            [i,j] = find([null, -speye(s,s); K, K]);
            this.JSparsityPattern = sparse(i,j,ones(size(i)),2*s,2*s);
        end
    end
end