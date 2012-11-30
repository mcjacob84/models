classdef MuscleFibreModel < models.BaseFullModel
% MuscleFibreModel: Model for a muscle fibre model composed of motoneuron, spindle and
% sarcomere array.
%
% @todo place citations
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
    end
    
    methods
        function this = MuscleFibreModel(N)
            if nargin < 1
                N = 100;
            end
            this.System = models.muscle.MuscleFibreSystem(this,N);
            
            this.T = .5; % [ms]
            this.dt = 1e-5; % [ms]
        end
    end
    
end