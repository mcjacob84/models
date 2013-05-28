classdef DynLinTimoshenkoModel < models.BaseFullModel
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
        
        % Flag that indicates if the nonlinear model is used.
        %
        % @type logical @default false 
        NonlinearModel = false;
        
        Supports;
        
        Loads;
        
        Materials;
        
        % The configuration file used for this model
        %
        % @type char
        ConfigFile;
        
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
    
    properties(SetAccess=private)
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
    end
    
    methods
        function this = DynLinTimoshenkoModel(cfgfile, nonlinear)
            % Creates a new timoshenko model
            %
            % Parameters:
            % cfgfile: The configuration file for the beams. @type char @default 'Simpel1.txt'
            % nonlinear: Flag to indicate if the linear or nonlinear model should be used.
            % @type logical @default false
            %
            % Return values:
            % this: The new model instance
            this = this@models.BaseFullModel;
            
            if nargin < 2
                this.NonlinearModel = false;
                if nargin < 1
                    %error('no config file given');
                    cfgfile = 'Simpel1.txt';
                    fprintf('No config file specified, using ''%s''\n',cfgfile);
                end
            else
                this.NonlinearModel = nonlinear;
            end
            this.Name = 'DynLin Timoschenko Beam';
            
            %% Load geometry from config file
            path = fileparts(mfilename('fullpath'));
            this.ConfigFile = fullfile(path,cfgfile);
            [this.Points, this.RO_raw, this.KR_raw, this.FH_raw, raw_mat, this.Supports, this.Loads] = this.read_file(this.ConfigFile);
            
            this.Materials = models.beam.Material.empty;
            for midx = 1:size(raw_mat,1)
                this.Materials(midx) = models.beam.Material(raw_mat(midx,:));
            end
            
            this.split_RO;
            this.split_KR;
            % Important: call before creating the system!
            this.preprocess_data;
            
            %% Model specifics
            this.T = 5;
            this.dt = .05;
            this.Data.useFileTrajectoryData;
            
            %% Internal setup
            this.System = models.beam.DynLinTimoshenkoSystem(this);
            
            % Subspace reduction
            s = spacereduction.PODGreedy;
            s.Eps = 1e-7;
            this.SpaceReducer = s;
            
            % No approximation - core fun is linear
            this.Approx = [];
            
            % No Parameters
            this.Sampler = sampling.GridSampler;
            %this.Sampler.Samples = 20;
            
            %% Function & ODE solver setup
            if this.NonlinearModel
                this.System.f = models.beam.DLTNonlinearCoreFun(this.System);
                
                % ODE Solver -> Use Matlab ode15i
%                 o = solvers.MLode15i;
%                 o.RelTol = 1e-3;
%                 o.AbsTol = 1e-3;
                
                o = solvers.FullyImplEuler(this);
                
                %o = solvers.MLWrapper(@ode45);
                %o.MaxStep = this.dt;
                
                this.ODESolver = o;
            else
                this.System.f = models.beam.DLTLinearCoreFun(this.System);
                this.ODESolver = solvers.LinearImplEuler(this);
            end
            
            % Train with all inputs
            this.TrainingInputs = 1:this.System.InputCount;
            
            % Setup JKerMor export
            je = export.JKerMorExport;
            je.JavaExportPackage = 'models.beam.dynlintimo';
            je.Short = this.Name;
            je.GeometryExportCallback = @this.exportGeometry;
            je.JaRMoSBaseSource = 'C:\Users\CreaByte\Documents\Uni\Software\JaRMoSBase\src';
            this.JKerMorExport = je;
        end
        
        function m = get.ColorMap(this)
            % col = colormap('hot');
            % Eigene Colormap erstellen (num_col Farben)
            % blau -> gr�n -> rot 
            % num_col = 128;
            % dc = [0:2/num_col:1]';
            % col = [0*dc dc 1-dc; dc 1-dc 0*dc];
            % blau -> gr�n -> gelb -> rot
            dc = [0:3/this.NumColors:1]';
            m = [0*dc dc 1-dc; dc 0*dc+1 0*dc; 0*dc+1 1-dc 0*dc];
        end
        
        function exportGeometry(this, f, folder, export)
            % Exports the model geometry to JKerMor
            %
            % To be done
            %
            % Parameters:
            % f: A file handle to write the model xml @type handle
            % folder: The target folder @type char
            % export: The exporter instance. @type export.JaRMoSExport
            
            fprintf(f,'\t<dimension>3</dimension>\n');
            fprintf(f,'\t<nodes>%d</nodes>\n',size(this.Points,1));
%             fprintf(f,'\t<fieldmapping>VERTEX</fieldmapping>\n');
            fprintf(f,'\t<hasFaces>false</hasFaces>\n');
            
            %% Nodes (including dirichlet points)
            usedpts = [];
            for i=1:length(this.Elements)
                el = this.Elements{i};
                usedpts = [usedpts el.PointsIdx];%#ok
            end
            usedpts = unique(usedpts);
            
            vert = single(reshape(this.Points(usedpts,:)',1,[]));
            export.saveRealVector(vert, 'vertices.bin', folder);
            
            % Dirichlet nodes
            fprintf(f,'\t<hasDirichletNodes>true</hasDirichletNodes>\n');
            export.saveRealVector(int16(this.dir_u), 'dir_nodes.bin', folder);
            export.saveRealVector(this.u_dir, 'dir_values.bin', folder);
            
            n = length(this.Elements);
            edges = int16(2*n);
            for i=1:n
                edges((2*i-1):(2*i)) = this.Elements{i}.PointsIdx;
            end
            export.saveRealVector(edges, 'edges.bin', folder);
        end
        
        plot(model, t, u);

        plotSingle(model, t, u, h);
    end
    
    %% Getter & Setters
    methods
        function set.BeamRefinementFactor(this, value)
            [this.Points, this.RO_raw, this.KR_raw] = this.read_file(this.ConfigFile);
            if value ~= this.RO_factor_global
                this.ModelData.SimCache.clearTrajectories;
            end
            this.RO_factor_global = value;
            this.split_RO;
            this.split_KR;
            % Update model
            this.preprocess_data;
            % Update system
            this.System.buildElementDependentComponents;
            % Update core function
            this.System.f.initialize;
        end
        
        function value = get.BeamRefinementFactor(this)
            value = this.RO_factor_global;
        end

        function set.CurvedBeamRefinementFactor(this, value)
            [this.Points, this.RO_raw, this.KR_raw] = this.read_file(this.ConfigFile);
            if value ~= this.KR_factor_global
                this.ModelData.SimCache.clearTrajectories;
            end
            this.KR_factor_global = value;
            this.split_RO;
            this.split_KR;
            % Update model
            this.preprocess_data;
            % Update system
            this.System.buildElementDependentComponents;
            % Update core function
            this.System.f.initialize;
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
    end
    
    methods(Static)
        function experiment1
            m = models.beam.DynLinTimoshenkoModel('Rohrleitungen z.txt');
            m.BeamRefinementFactor = 5;
            m.TrainingInputs = 1:m.System.InputCount;
            m.T = 10;
            m.dt = .1;
            m.System.Params(1).Range = [0 2];
            m.System.Params(1).Desired = 10;
            m.System.Params(2).Range = [0 .1];
            m.System.Params(2).Desired = 10;
            m.System.Params(3).Range = 20*[-1 1];
            m.System.Params(3).Desired = 10;
            m.offlineGenerations;
            r = m.buildReducedModel;
            save experiment1 m r;
        end
    end
    
end