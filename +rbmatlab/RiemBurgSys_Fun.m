classdef RiemBurgSys_Fun < models.BaseDynSystem & dscomponents.ACoreFun
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
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods
        function this = RiemBurgSys_Fun(model)
            this = this@models.BaseDynSystem(model);
            this = this@dscomponents.ACoreFun;
            
            this.x0 = dscomponents.PointerInitialValue(@(mu)this.getx0(mu));
            
            % DS-Components
            this.f = this;
            this.B = [];
            this.Inputs = {};
            
            % Parameters
            this.addParam('ULeft', [.3, .3], 1);
            this.addParam('URight', [.1, 1], 30);
            this.addParam('xFlux', [-.5, -.5], 1);
            
            % Output conversion
            %c = zeros(1, params.xnumintervals*params.ynumintervals);
            %c(round(.75*params.xnumintervals)) = 1;
            %this.C = dscomponents.PointerOutputConv(@(t,mu)c,false);
        end
        
        function y = evaluateCoreFun(this, x, t, mu)
            rbm = this.Model.RBMModel;
            rbm = rbm.set_mu(rbm,mu);
            rbm = rbm.set_time(rbm,t);
            y = -rbm.L_E_local_ptr(rbm, this.Model.RBMDataCont.RBMData, x, []);
        end
        
        function x0 = getx0(this, mu)
            rbm = this.Model.RBMModel;
            rbm = rbm.set_mu(rbm,mu);
            % initial values by midpoint evaluation
            x0 = rbm.init_values_algorithm(rbm,this.Model.RBMDataCont.RBMData);
        end
    end
    
end