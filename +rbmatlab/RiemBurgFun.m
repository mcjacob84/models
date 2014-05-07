classdef RiemBurgFun < dscomponents.ACoreFun
% RiemBurgSys: 
%
%
%
% @author Daniel Wirtz @date 2011-11-15
%
% @new{0,6,dw,2011-11-15} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    methods
        function this = RiemBurgFun(sys)
            this = this@dscomponents.ACoreFun(sys);
            this.TimeDependent = true;
        end
        
        function y = evaluateCoreFun(this, x, t, mu)
            rbm = this.Model.RBMModel;
            rbm = rbm.set_mu(rbm,mu);
            rbm = rbm.set_time(rbm,t);
            y = -rbm.L_E_local_ptr(rbm, this.Model.RBMDataCont.RBMData, x, []);
        end
    end
    
end