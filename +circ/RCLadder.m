classdef RCLadder < models.BaseFullModel
% RCLadder: 
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
    
    properties
        % The target dimension of the RC ladder circuit
        %
        % @propclass{data}
        Dims;
    end
    
    methods
        
        function this = RCLadder(dims)
            if nargin == 0
                dims = 30;
            end
            this.Dims = dims;
            
            this.registerProps('Dims');
            
            this.Name = 'RC Ladder circuit model';
            
            this.T = 1;
            this.dt = .005;
            this.tau = 1;
            
            this.System = models.circ.RCLadderSys(this);
            
            this.Sampler = [];
            
            a = approx.AdaptiveCompWiseKernelApprox;
            a.TimeKernel = kernels.NoKernel;
            a.ParamKernel = kernels.NoKernel;
            a.SystemKernel = kernels.GaussKernel(.5);
            a.MaxRelErr = 1e-5;
            a.MaxAbsErrFactor = 1e-3;
            a.NumGammas = 10;
            this.Approx = a;
            
            s = spacereduction.PODReducer;
            s.Mode = 'abs';
            s.Value = 3;
            s.UseSVDS = false;
            this.SpaceReducer = s;
            
            s = solvers.ode.ExplEuler;
            s.MaxStep = [];
            this.ODESolver = s;
        end
    end
    
end