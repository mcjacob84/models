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

    properties(Constant)
        % No heat modeling at the current stage
        withHeat = false;
    end
    
    properties(SetObservable)        
        % Plot enhancement factor
        %
        % @propclass{optional} Scales the displacement
        %
        % @type double @default 1
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
        
        % Flag that determines if the nonlinear model should be used.
        %
        % @propclass{optional} Switches between model types
        %
        % @default false @type logical
        NonlinearModel;
    end
    
    properties(SetAccess=private)
        % 3D coordinates of the beam end points
        %
        % @type matrix
        Points;
        
        % Cell array of all model structure elements.
        % 
        % Fields:
        % p: The index of the start and end point of the beam, defined in
        % property Points.
        %
        % @type struct
        Elements;
        
        Supports;
        
        Loads;
        
        Materials;
        
        data;
        
        % Vector for neumann boundary conditions
        f_neum;
        
        dir_u;
        dir_T;
        
        % The indices in the global state space vector of all points
        % including dirichlet points (dim=7: 3location, 3velocity & heat) 
        free;
        
        % Extracted Dirichlet values full u vector (using dir_u)
        u_dir;
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
        
        RO_factor_global = 1;
        
        KR_factor_global = 1;
        
        fNonlin = false;
    end
    
    methods
        function this = DynLinTimoshenkoModel(cfgfile)
            if nargin == 0
                %error('no config file given');
                cfgfile = 'Simpel1.txt';
                fprintf('No config file specified, using ''Simpel1.txt''\n');
            end
            
            this = this@models.BaseFullModel;
            this.Name = 'DynLin Timoschenko Beam';
            this.JavaExportPackage = 'models.beam.dynlintimo';
            
            %% Load geometry from config file
            path = fileparts(mfilename('fullpath'));
            cfile = fullfile(path,cfgfile);
            [this.Points, this.RO_raw, this.KR_raw, this.FH_raw, raw_mat, this.Supports, this.Loads] = this.read_file(cfile);
            
            this.Materials = models.beam.Material.empty;
            for midx = 1:size(raw_mat,1)
                this.Materials(midx) = models.beam.Material(raw_mat(midx,:));
            end
            
            this.split_RO;
            this.split_KR;
            this.preprocess_data;
            
            %% Model specifics
            this.T = 5;
            this.dt = .05;
            
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
            
            %% Call update of model type (sets the system's f function and chooses an ODE solver)
            this.updateModelType(this.NonlinearModel);
            
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
        
        function v = get.NonlinearModel(this)
            v = this.fNonlin;
        end
        
        function set.NonlinearModel(this, v)
            if ~islogical(v) || ~isscalar(v)
                error('NonlinearModel must be a scalar logical');
            end
            if (this.fNonlin ~= v)
               this.updateModelType(v);
            end
            this.fNonlin = v;
        end
    end
    
    methods(Access=private)
        function updateModelType(this, nonlin)
            % Updates the model type (switch for nonlinear/linear)
            %
            % Automatically chooses the correct CoreFun and a suitable ODE
            % solver.
            s = this.System;
            if nonlin
                s.f = models.beam.DLTNonlinearCoreFun(s);
                % ODE Solver -> Use Matlab ode15i
                o = solvers.ode.MLode15i;
                o.RelTol = 1e-3;
                o.AbsTol = 1e-3;
                
                %o = solvers.ode.MLWrapper(@ode45);
                %o.MaxStep = this.dt;
                
                this.ODESolver = o;
            else
                s.f = models.beam.DLTLinearCoreFun(s);
                this.ODESolver = solvers.ode.LinearImplEuler(this);
            end
        end
        
        % Reads the config file
        [Points, RO_raw, KR, FH, mat, lager, lasten] = read_file(this, file);
        
        [data, RO, KR, FH, p] = preprocess_data(p, RO_raw, KR_raw, FH_raw, c, supports, loads, gravity, RO_factor_global, KR_factor_global);
        
        split_RO(this);
        
        split_KR(this);
                
        plot_circle(this, N, T, T1, T2, T_Fren, R, angle, B, pc, u1, u2, col1, col2, plot_options);
    end
    
end