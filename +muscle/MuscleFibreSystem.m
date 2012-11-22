classdef MuscleFibreSystem < models.BaseDynSystem
% MuscleFibreSystem: The global dynamical system used within the MuscleFibreModel
%
% Contains the FibreDynamics, Linear Diffusion and InputConv components as well as a
% ConstInitialValue.
%
% @author Daniel Wirtz @date 2012-11-22
%
% @new{0,7,dw,2012-11-22} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        % The number of sarcomer cells in this fibre model
        N;
    end
    
    methods
        function this = MuscleFibreSystem(model, N)
            % Call superclass constructor
            this = this@models.BaseDynSystem(model);
            % Set local variables
            this.N = N;
            
            %% Set system components
            % Core nonlinearity
            this.f = models.muscle.FibreDynamics(N);
            
            % Linear diffusion part
            this.A = dscomponents.LinearCoreFun(this.assembleA);
            
            % Linear input B for motoneuron
            % TODO B (dscomponents.LinearInputConv)
            
            % Constant initial values
            this.x0 = dscomponents.ConstInitialValue(this.initStates);
        end
    end
    
    methods(Access=private)
        function A = assembleA(this)
               % this.N 
               A = 0;
        end
        
        function x0 = initStates(this)
             x0 = 0;
             return;
             x0 = [this.initNeuroStates;...
                this.initSpindleStates;...
                this.initSarcomerStates(this.N)];
        end
        
        function STATES = initNeuroStates(~)
            STATES = zeros(6,1);
        end
    end
    
end