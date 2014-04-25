classdef BaseCoreFun < dscomponents.ACompEvalCoreFun
% BaseCoreFun: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2013-03-14
%
% @new{0,7,dw,2013-03-14} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties(Constant)
        % The time in seconds needed before the activation function reaches its maximum value
        % and decays back to zero.
        %
        % Only applies for ActivationFunType == 2
        %
        % @default 30 @type double
        %
        % See also: MaxActivationTime ActivationFunType
        ActivationTransitionTime = 10; %[s]
        
        % The maximum time in seconds that spans the support of the piecewise activation
        % function.
        % It is composed of gaussian-shaped increase/decrease ends and a constant one in
        % between. For the default settings, we would have a maximum duration of
        % 500s-2*30s=440s for the activation rate of level one.
        %
        % The model parameter "atime" determines the level-one-activation time linearly between
        % [2*ActivationTransitionTime, MaxActivationTime]
        %
        % The activation function always starts at t=0.
        %
        % Only applies for ActivationFunType == 2
        %
        % @default 500 @type double
        %
        % See also: ActivationTransitionTime ActivationFunType
        MaxActivationTime = 400; %[s]
    end

    properties(Dependent)
        % Type of the activation function.
        %
        % Included for backward compatibility with older simulations.
        %
        % Admissible values:
        % 1: First version with gaussian activation over 300s
        % 2: Parameterized version
        %
        % @default 1 @type integer
        ActivationFunType;
    end
    
    properties(Access=private)
        gaussian;
        fAFT = 1;
    end
    
    methods
        function this = BaseCoreFun(dynsys)
            this = this@dscomponents.ACompEvalCoreFun(dynsys);
            this.TimeDependent = true;
            
            this.ActivationFunType = 1;
        end
        
        function copy = clone(this, copy)
            % Call superclass method
            copy = clone@dscomponents.ACompEvalCoreFun(this, copy);
            % Sets gaussian etc
            copy.ActivationFunType = this.ActivationFunType;
        end
        
        function evaluateCoreFun(varargin)
            error('dont call me (direct overload of evaluate for efficiency');
        end
        
        function plotActivationFun(this, mu, pm)
            if nargin < 3
                pm = PlotManager;
                pm.LeaveOpen = true;
                if nargin < 2
                    mu = this.System.Model.getRandomParam;
                    
                end
            end
            
            h = pm.nextPlot('activation_fun','External input activation function','time','factor');
            t = linspace(0,min(this.System.Model.T,this.MaxActivationTime*1.2),2000);
            plot(h,t,this.activationFun(t/this.System.Model.tau,mu));
            
            if nargin < 3
                pm.done;
            end
        end
        
        function set.ActivationFunType(this, value)
            % Activation function setup
            if value == 1
                % \gamma=28 is such that we have K(0,27)~<.4 (27=150/tau)
                k = kernels.GaussKernel(28.206364723698);
            else
                k = kernels.GaussKernel;
                k.setGammaForDistance(this.ActivationTransitionTime/this.System.Model.tau,.001);
            end
            this.gaussian = k;
            this.fAFT = value;
        end
        
        function value = get.ActivationFunType(this)
            value = this.fAFT;
        end
    end
    
    methods(Access=protected)
        function f = activationFun(this, t, mu)
            g = this.gaussian;
            if this.ActivationFunType == 1
                f = (g.evaluateScalar(t-27)-.4).*(t<=54)/.6;
            else
                tau = this.System.Model.tau;
                ts = this.ActivationTransitionTime/tau;
                te = ts+(mu(3,:)*(this.MaxActivationTime-2*this.ActivationTransitionTime))/tau;
                f = (g.evaluateScalar(t-ts)-.001).*(t<=ts)/.999 ...
                    + (t>ts).*(t<=te) ...
                    + (g.evaluateScalar(t-te)-.001).*(t>te).*(t<=te+ts)/.999;
            end
        end
    end
    
    methods(Static, Access=protected)
        function obj = loadobj(obj, from)
            if nargin == 2
                if isfield(from,'System')
                    obj.System = from.System;
                end
                % Sets gaussian etc
                if isfield(from,'fAFT')
                    AFT = from.fAFT;
                else
                    AFT = 1;
                end
                obj.ActivationFunType = AFT;
                obj = loadobj@dscomponents.ACompEvalCoreFun(obj, from);
            else
                obj = loadobj@dscomponents.ACompEvalCoreFun(obj);
            end
        end
    end
    
end