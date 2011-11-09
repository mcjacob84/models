classdef DynLinTimoshenkoCoreFun < dscomponents.LinearCoreFun
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
    
    properties
    end
    
    methods
        function this = DynLinTimoshenkoCoreFun(model)
            dim = model.dim;
            A = (rand(dim,dim)-.7) * .01;
            this = this@dscomponents.LinearCoreFun(A);
        end
        
        function prepareConstants(mu, inputidx)
            % do constant preps, i.e. assign A
        end
    end
    
end