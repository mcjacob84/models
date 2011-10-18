classdef DynLinTimoshenkoModel < models.BaseFullModel & export.JKerMorExportable
% DynLinTimoshenkoModel: 
%
%
%
% @author Daniel Wirtz @date 2011-09-20
%
% @new{0,5,dw,2011-09-20} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetObservable)        
        % Plot enhancement factor
        %
        % @propclass{optional} Scales the displacement
        %
        PlotFactor = 20;
        
        % The beam line width for plotting
        %
        % @propclass{optional} Visualization setting
        BeamLineWidth = 3;
        
        % Number of colors for plotting
        %
        % @propclass{optional}
        %
        % @see ColorMap
        NumColors = 128;
    end
    
    properties(SetAccess=private)
        % The system's full dimension
        dim;
        
        % 3D coordinates of the beam end points
        Points;
        
        % Struct for each beam with several fields.
        % 
        % Fields:
        % p: The index of the start and end point of the beam, defined in
        % property Points.
        Beams;
    end
    
    properties(Dependent, SetAccess=private)
        % The color map to use. Value is hardcoded so far.
        ColorMap;
    end
    
    properties(Access=private)
        % Maximum temperature for plotting
        maxTemp = 0;
        
        % Maximum temperature for plotting
        minTemp = 0;
        
        % ColorBar handle
        cbh;
    end
    
    methods
        function this = DynLinTimoshenkoModel
            this = this@models.BaseFullModel;
            this.Name = 'DynLin Timoschenko Beam';
            this.JavaExportPackage = 'models.beam.dynlintimo';
            
            %% Load geometry from config file
            path = fileparts(mfilename('fullpath'));
            cfile = fullfile(path,'Config.txt');
            [this.Points, RO_raw, KR, FH, mat, lager, lasten] = this.read_file(cfile);%#ok
            nBeams = size(RO_raw, 1);
            this.Beams = struct('p',{});
            for i = 1:nBeams
                this.Beams(i).p = [RO_raw(i,2) RO_raw(i,3)];
            end
            
            %% Model specifics
            this.T = 5;
            this.dt = .05;
            % Dimensions: 1..3 spatial (x,y,z), 4..6 velocity (x,y,z) 7
            % temperature
            this.dim = 7*size(this.Points,1);
            
            %% Internal setup
            this.System = models.beam.DynLinTimoshenkoSystem(this);
            
            % Subspace reduction
            s = spacereduction.TrajectoryGreedy;
            s.Eps = 1e-7;
            this.SpaceReducer = s;
            
            % No approximation - core fun is linear
            this.Approx = [];
            
            % No Parameters
            this.Sampler = [];%sampling.RandomSampler;
            %this.Sampler.Samples = 20;
            
            % ODE Solver -> Use Matlab ode15i
            o = solvers.ode.MLode15i;
            o.RelTol = 1e-4;
            o.AbsTol = 1e-5;
            o.MaxStep = [];
            this.ODESolver = o;
            
            % Train with all inputs
            this.TrainingInputs = [1 2 3];
        end
        
        function m = get.ColorMap(this)
            % col = colormap('hot');
            % Eigene Colormap erstellen (num_col Farben)
            % blau -> grün -> rot 
            % num_col = 128;
            % dc = [0:2/num_col:1]';
            % col = [0*dc dc 1-dc; dc 1-dc 0*dc];
            % blau -> grün -> gelb -> rot
            dc = [0:3/this.NumColors:1]';
            m = [0*dc dc 1-dc; dc 0*dc+1 0*dc; 0*dc+1 1-dc 0*dc];
        end
        
        function exportGeometry(this, f)
            % To be done
            
        end
    end
    
    methods(Access=private)
        [Points, RO_raw, KR, FH, mat, lager, lasten] = read_file(this, file);
    end
    
end