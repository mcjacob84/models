classdef DLTLinearCoreFun < models.beam.DLTBaseCoreFun & dscomponents.AffLinCoreFun
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
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    methods
        function this = DLTLinearCoreFun(system)
            this = this@models.beam.DLTBaseCoreFun(system);
            this.initialize;
            % Set class for export
            this.CoeffClass = 'CoreFunCoeffs';
        end
        
        function initialize(this)
            initialize@models.beam.DLTBaseCoreFun(this);
            
            m = this.sys.Model;
            
            
            
            
        end
    end
end