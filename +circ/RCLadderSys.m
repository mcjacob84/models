classdef RCLadderSys < models.BaseDynSystem
% RCLadderSys: 
%
%
%
% @author Daniel Wirtz @date 2011-04-29
%
% @new{0,3,dw,2011-04-29} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods
        function this = RCLadderSys(model)
            this = this@models.BaseDynSystem(model);
            
            dim = this.Model.Dims;
            
            B = zeros(dim,1);
            B(1) = 1;
            
            this.B = dscomponents.LinearInputConv(B);
            this.C = dscomponents.LinearOutputConv(B');
            
            this.f = models.circ.RCLadderFun;
            
            %this.Inputs = {@(t)exp(-t),@(t)sin(t),@(t)0};
            this.Inputs = {@(t)exp(-t),@(t)(cos(2*pi*t/10) + 1)/2};
            
            this.MaxTimestep = [];
            this.StateScaling = 1;            
        end
        
        function x = x0(this, mu)%#ok
            x = zeros(this.Model.Dims,1);
        end
    end
    
    methods(Access=protected)
        function validateModel(this, model)%#ok
            % Validates if the model to be set is a valid BaseModel at
            % least.
            % Extracting this function out of the setter enables subclasses
            % to further restrict the objects that may be passed, as is
            % being done in models.ReducedSystem, for example.
            if ~isa(model, 'models.circ.RCLadder')
                error('The Model property has to be a child of models.circ.RCLadder');
            end
        end
    end
    
end