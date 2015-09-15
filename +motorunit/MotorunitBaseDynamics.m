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
        
        function f = getLinkFactor(this, y)
            % scalar link factor motoneuron to middle sarcomere
            %
            % See also: FibreDynamics.evaluate
            % Parameters:
            % y: @type rowvec<double>  second state from motoneuron
            % submodel
            f = this.MSLink_MinFactor + exp(-(y-this.MSLink_MaxFactorSignal).^2/150)...
                *(this.MSLink_MaxFactor-this.MSLink_MinFactor);
        end
    end
    
end