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
            
            this.x0 = dscomponents.ConstInitialValue(zeros(this.Model.Dims,1));
            
            B = zeros(dim,1);
            B(1) = 1;
            this.B = dscomponents.LinearInputConv(B);
                
            e = ones(dim,1);
            A = spdiags([e -2*e e],-1:1,dim,dim);
            A(1,:) = -A(1,:);
            A(end,end) = -1;
            this.A = dscomponents.LinearCoreFun(A);
            
            this.C = dscomponents.LinearOutputConv(B');
            
            this.f = models.circ.RCLadderFun(dim);
            
            this.Inputs = {@(t)exp(-t),@(t)(cos(2*pi*t/10) + 1)/2,...
                @(t)double((t-.3)>0),@(t)(cos(pi*t)+1)/(t+1),...
                @(t)(cos(pi*t*.5)+1)/(t+1),@(t)2*(t>1.3),@(t)(cos(3*pi*t)+1)/(2*t+1)};
            
            this.MaxTimestep = [];
            this.StateScaling = 1;            
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