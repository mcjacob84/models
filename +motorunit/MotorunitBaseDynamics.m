classdef MotorunitBaseDynamics < dscomponents.ACoreFun
% MotorunitBaseDynamics: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2015-09-15
%
% @new{0,7,dw,2015-09-15} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/
% - \c Documentation http://www.morepas.org/software/kermor/
% - \c License @ref licensing
    
    properties
        % The `V_s` value of the motoneuron input at which the MSLink_MaxFactor should be attained
        %
        % @type double @default 40
        MSLink_MaxFactorSignal = 40;
        
        % The maximal factor with which the `V_s` value of the motoneuron should be amplified
        % when added to the sarcomere equations
        %
        % @type double @default 7
        MSLink_MaxFactor = 7;
        
        % The minimal factor at which the `V_s` value of the motoneuron should be amplified
        % when added to the sarcomere equations.
        %
        % Between both limits, an exponential Gaussian weight is applied with a certain radius,
        % so that the amplification factor has its minimum value at zero and maximum value at
        % MSLink_MaxFactorSignal.
        %
        % See also FibreDynamics.getLinkFactor
        %
        % @type double @default .3
        MSLink_MinFactor = .3;
    end
    
    methods
        function this = MotorunitBaseDynamics(sys)
            this = this@dscomponents.ACoreFun(sys);
        end
        
        function factor = getLinkFactor(this, signal)
            % scalar link factor motoneuron to middle sarcomere
            %
            % See also: Dynamics.evaluate
            %
            % Parameters:
            % signal: @type rowvec<double>  second state from motoneuron
            % submodel            
            factor = this.MSLink_MaxFactor*ones(1,length(signal));
            dynfac = signal < this.MSLink_MaxFactorSignal;
            factor(dynfac) = this.MSLink_MinFactor + exp(-(signal(dynfac)-this.MSLink_MaxFactorSignal).^2/150)...
                *(this.MSLink_MaxFactor-this.MSLink_MinFactor);
        end
        
        function pm = plotMotoSacroLinkFactorCurve(this)
            x = 0:.1:80;
            pm = PlotManager;
            pm.LeaveOpen = true;
            h = pm.nextPlot('moto_sarco_link_factor','Factor for motoneuro to sarcomere link','Moto V_s','Factor');
            fx = this.MSLink_MaxFactor*ones(1,length(x));
            dynfac = x < this.MSLink_MaxFactorSignal;
            fx(dynfac) = this.getLinkFactor(x(dynfac));
            plot(h,x,fx);
            pm.done;
        end
    end
    
end