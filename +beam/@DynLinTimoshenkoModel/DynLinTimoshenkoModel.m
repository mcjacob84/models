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
        % @type double @default 20
        PlotFactor = 20;
        
        % The beam line width for plotting
        %
        % @propclass{optional} Visualization setting
        %
        % @type double @default 3
        BeamLineWidth = 3;
        
        % Number of colors for plotting
        %
        % @propclass{optional}
        %
        % @type integer @default 128
        %
        % @see ColorMap
        NumColors = 128;
    end
    
    properties(SetObservable, Dependent)
        BeamRefinementFactor;
        
        CurvedBeamRefinementFactor;
    end
    
    properties(SetAccess=private)
        % The system's full dimension
        %
        % @type integer
        dim;
        
        % 3D coordinates of the beam end points
        %
        % @type matrix
        Points;
        
        % Struct for each beam with several fields.
        % 
        % Fields:
        % p: The index of the start and end point of the beam, defined in
        % property Points.
        %
        % @type struct
        Beams;
        
        Supports;
        
        Loads;
        
        RO;
        KR;
        FH;
        
        data;
    end
    
    properties(Dependent, SetAccess=private)
        % The color map to use. Value is hardcoded so far.
        %
        % @type matrix
        ColorMap;
    end
    
    properties(Access=private)
        % Maximum temperature for plotting
        maxTemp = 0;
        
        % Maximum temperature for plotting
        minTemp = 0;
        
        % ColorBar handle
        cbh;
        
        RO_raw;
        
        KR_raw;
        
        FH_raw;
        
        RO_factor_global = 5;
        
        KR_factor_global = 5;
        
        mat;
    end
    
    methods
        function this = DynLinTimoshenkoModel(cfgfile)
            if nargin == 0
                error('no config file given');
            end
            
            this = this@models.BaseFullModel;
            this.Name = 'DynLin Timoschenko Beam';
            this.JavaExportPackage = 'models.beam.dynlintimo';
            
            %% Load geometry from config file
            path = fileparts(mfilename('fullpath'));
            cfile = fullfile(path,cfgfile);
            [this.Points, this.RO_raw, this.KR_raw, this.FH_raw, this.mat, this.Supports, this.Loads] = this.read_file(cfile);
            this.split_RO;
            this.split_KR;
            this.preprocess_data;
            
            %% Model specifics
            this.T = 5;
            this.dt = .05;
            % Dimensions: 1..3 spatial (x,y,z), 4..6 velocity (x,y,z) 7
            % temperature
            this.dim = 14*this.data.num_knots;
            
            %% Internal setup
            this.System = models.beam.DynLinTimoshenkoSystem(this);
            
            % Subspace reduction
            s = spacereduction.PODGreedy;
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
%             this.TrainingInputs = [1 2 3];
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
            % Exports the model geometry to JKerMor
            %
            % To be done
            %
            % Parameters:
            % f: A file handle to write the model xml @type handle
        end
        
        plot(model, t, u);

        plotSingle(model, t, u);
        
        N = circle_shape_functions(this, R, s, B);
    end
    
    %% Getter & Setters
    methods
        function set.BeamRefinementFactor(this, value)
            this.RO_factor_global = value;
            this.split_RO;
            this.preprocess_data;
        end
        
        function value = get.BeamRefinementFactor(this)
            value = this.RO_factor_global;
        end

        % Todo: für KR!
        function set.CurvedBeamRefinementFactor(this, value)
            this.KR_factor_global = value;
            this.split_KR;
            this.preprocess_data;
        end
        
        function value = get.CurvedBeamRefinementFactor(this)
            value = this.KR_factor_global;
        end
    end
    
    methods(Access=private)
        % Reads the config file
        [Points, RO_raw, KR, FH, mat, lager, lasten] = read_file(this, file);
        
        [data, RO, KR, FH, p] = preprocess_data(p, RO_raw, KR_raw, FH_raw, c, supports, loads, gravity, RO_factor_global, KR_factor_global);
        
        split_RO(this);
        
        split_KR(this);
        
        plot_beam(this, split, T, c, p1, p2, u1, u2, col1, col2, plot_options);
        
        plot_circle(this, N, T, T1, T2, T_Fren, R, angle, B, pc, u1, u2, col1, col2, plot_options);
        
        N = beam_shape_functions_derivative(this, s, L, c);
        
        N = beam_shape_functions(this, s, L, c);
        
        B = circle_connect_matrix(this, R, L);
    end
    
end