classdef PCDModel < models.BaseFullModel
% Base class for both 1D and 2D pcd models
%
% @author Daniel Wirtz @date 2011-03-16
%
% @change{0,6,dw,2011-11-27} Moved the static reduction experiments to their respective xD
% systems.
%
% @new{0,5,dw,2011-11-02} Added many static reduction experiment cases to the model.
%
% @new{0,3,dw,2011-03-16} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(Constant)
        % Typical diffusion rate for Caspase-8
        d1 = 1.8e-11; %[m^2/s]
        
        % Typical diffusion rate for Caspase-3
        d2 = 1.86e-11; %[m^2/s]
        
        % Typical diffusion rate for Pro-Caspase-8
        d3 = 1.89e-11; %[m^2/s]
        
        % Typical diffusion rate for Pro-Caspase-3
        d4 = 2.27e-11; %[m^2/s]
        
        % Typical cell length (from 1D)
        L = 1e-5; %[m]
    end
    
    properties(Dependent)
        Dimension;
    end
    
    methods
        function this = PCDModel(dim)
            % Creates a new instance of the PCDModel
            %
            % Parameters:
            % dim: The dimension to use @default 1
            
            if nargin == 0
                dim = 1;
            end
            
            this.T = 6; %[s]
            this.dt = .01; %[s]
            % time scaling
            this.tau = this.L^2/this.d1;
            
            this.Data.TrajectoryData = data.FileTrajectoryData(this.Data);
            
%             s = sampling.RandomSampler;
%             s.Samples = 10;
            s = sampling.GridSampler;
            s.Spacing = 'log';
            this.Sampler = s;
            
%             this.ODESolver = solvers.MLWrapper(@ode23);
%             this.ODESolver = solvers.ExplEuler(this.dt);
            o = solvers.MLode15i;
            o.AbsTol = 1e-6;
            o.RelTol = 1e-5;
            o.MaxStep = [];
            this.ODESolver = o;
            
            switch dim
                case 3
                    s = models.pcd.PCDSystem3D(this);
                    this.Name = 'Programmed Cell Death 3D';
                case 2 
                    s = models.pcd.PCDSystem2D(this);
                    this.Name = 'Programmed Cell Death 2D';
                otherwise
                    s = models.pcd.PCDSystem1D(this);
                    this.Name = 'Programmed Cell Death 1D';
            end
            this.System = s;
            
            % Space reduction setup
            sr = spacereduction.PODGreedy;
            sr.Eps = 1e-5;
            this.SpaceReducer = sr;
            
            this.Approx = [];
        end
        
        function plot(this, varargin)
            % Overrides standard method and forwards to the system's plot
            % function. (they are 1D and 2D)
            this.System.plot(this, varargin{:});
        end
        
        function plotState(this, varargin)
            % Overrides standard method and forwards to the system's plot
            % function. (they are 1D and 2D)
            this.System.plotState(this, varargin{:});
        end
        
        function value = get.Dimension(this)
            value = this.System.Dims;
        end
    end
end

