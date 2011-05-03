classdef RCLadder < models.BaseFullModel
% RCLadder: Model of a nonlinear resistor with independent current source
%
% This model has been used as benchmark in many papers. This implementation follows the details in
% [BS06] - Bai, Z., Skoogh, D.: A projection method for model reduction of bilinear dynamical systems.
% Linear Algebra and its Applications 415(2-3), 406?425 (2006)
%
% This model is also used in several papers dealing with MOR of nonlinear systems, i.e.
% [Re03] Rewienski, M.: A trajectory piecewise-linear approach to model order reduction of nonlinear
% dynamical systems. Ph.D. thesis, Citeseer (2003)
% or
% [CI04] M. Condon and R. Ivanov. Empirical balanced truncation of nonlinear systems. Journal of
% Nonlinear Science, 14:405?414, 2004. 10.1007/s00332-004-0617-5.
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
            this.dt = .0025;
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
            s.MaxStep = .005; % Stability constraint due to diffusion term
            this.ODESolver = s;
        end
    end
    
end