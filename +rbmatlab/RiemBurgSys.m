classdef RiemBurgSys < models.BaseFirstOrderSystem
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
        function this = RiemBurgSys(model)
            this = this@models.BaseFirstOrderSystem(model);
            this.TimeDependent = true;
            
            this.x0 = dscomponents.PointerInitialValue(@(mu)this.getx0(mu));
            
            % DS-Components
            this.f = models.rbmatlab.RiemBurgFun(this);
            this.B = [];
            this.Inputs = {};
            
            % Parameters
            this.addParam('ULeft', .3);
            this.addParam('URight', .5, 'Range', [.1, 1], 'Desired', 30);
            this.addParam('xFlux', -.5);
            
            % Output conversion
            %c = zeros(1, params.xnumintervals*params.ynumintervals);
            %c(round(.75*params.xnumintervals)) = 1;
            %this.C = dscomponents.PointerOutputConv(@(t,mu)c,false);
        end
        
        function x0 = getx0(this, mu)
            rbm = this.Model.RBMModel;
            rbm = rbm.set_mu(rbm,mu);
            % initial values by midpoint evaluation
            x0 = rbm.init_values_algorithm(rbm,this.Model.RBMDataCont.RBMData);
        end
    end
    
end