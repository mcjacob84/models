classdef Model < models.BaseFullModel
    % Model: Electromyography model
    %
    % Model representing the second equation in the bidomain equations.
    % `\nabla\cdot[(\sigma_i+\sigma_e)\nabla\phi_e] + \nabla\cdot[\sigma_i\nabla\phi_e] = 0`
    %
    % @author Timm Strecker @date 2014-03-24
    %
    % @new{0,7,ts,2014-03-24} Added this class.
    %
    properties
        % Seed that can be used by random number generator instances in order to enable result
        % reproduction.
        %
        % @type int @default 1
        RandSeed = 1;
    end
    
    properties(Dependent)   % properties related to the motor unit configuration. 
        
        % distribution of motor units over cross section
        %
        % @default ones(dim(2),zdim_m)
        MUdistribution;
        
        % the number of motor units
        %
        % @default 1
        numMU;
        
        % Flag indicating whether to use custom firing times of the motor units.
        % Otherwise, the firing times are generated automatically by using the parameter mu.
        FiringTimesMode;
        
        % Position (index) of the neuromuscular junction of each fibre. By
        % default randomly anywhere on the fibres.
        % matrix<integer> of dimension dim(2) \times zdim_m
        neuronpos;
        
        % Cell array of structs containing information about the motor
        % units. There are the fields
        %'type': specifying the type (slow-twitch, fast-twitch, etc.),
        % 'firing_times': containing the times when the motoneuron is firing,
        % 'centre' and 'radius': optional fields that define a circle
        % of radius radius around centre. All fibres of the MU are
        % distributed inside this circle.
        MUs;
        
        % the conductivity tensor in intracellular space
        %
        % @type matrix<double>
        sigmai = diag([8.93 0.893 0.893]); %eye(3);
        
        % the conductivity tensor in extracellular space
        %
        % @type matrix<double>
        sigmae = diag([6.7 6.7 6.7]);  %eye(3);
        
        % the conductivity tensor in the body (skin/fat)
        %
        % @type matrix<double>
        sigmao;
    end
    
    properties(SetAccess=private)
        OriginalV;
    end
    
    methods
        function this = Model(musclegeometry, dim)
            % Creates a new Electromyography model
            %
            % musclegeometry: vector with 4 entries representing length,
            % width, muscle height and height of the fat/skin layer (in cm). If
            % there are only 3 entries, skinheight is set to zero
            % Dimension of the discretization in all 3 dimensions.
            if nargin < 3
                if nargin < 2
                    dim = [40;20;13];
                    if nargin < 1
                        musclegeometry = 0.1*[40;20;10;2];
                    end
                end
            end
            if size(musclegeometry, 1) == 1
                musclegeometry = musclegeometry';
            end
            if size(dim, 1) == 1
                dim = dim';
            end
            if size(musclegeometry,1) == 3
                musclegeometry = [musclegeometry; 0];
            elseif size(musclegeometry,1) ~= 4
                error('Muscle geometry must have 3 or 4 entries.')
            end
            if size(dim, 1) ~= 3
                error('Dim must have 3 entries (being grid dimension in x-, y-, and z-direction).')
            end
            
            this.Name = sprintf('EMG model geometry = [%g,%g,%g,%g] discretization = [%g,%g,%g]',musclegeometry, dim);
            this.System = models.emg.System(this, musclegeometry, dim);
            %this.TrainingInputs = 1;
            this.SaveTag = sprintf('EMG_geom[%g,%g,%g,%g]_discr[%g,%g,%g]',musclegeometry, dim);
            this.Data = data.ModelData(this);
            this.isStatic = true;
            this.DefaultMu = 5;
            
            d = approx.DEIM(this.System);
            d.MaxOrder = 2000;
            this.Approx = d;
            s = spacereduction.PODReducer;
            s.Value=1e-5;
            s.Mode='eps';
            this.SpaceReducer=s;
            
            this.EnableTrajectoryCaching = true;
            
            this.T = 100;
            this.dt = 1;
        end
        
        function distributeMUs(this, varargin)
            % Distributes number many motor units over muscle-crosssection.
            % All fibres of each MU lie within a circle. By default, the
            % circles are distributed approximately uniformly over the
            % muscle-crosssection such that every fibre can be assigned to
            % at least one MU. The optional arguments centers and radii can
            % be used to define custom circles.
            %
            % Parameters:
            % number: The number of motor units
            % @type integer
            % centers: centers of the circles @type 2 \times number matrix<double>
            % radii: radii of the circles @type vector<double> of length
            % number or scalar
            this.System.f.distributeMUs(varargin{:});
        end
        
        function initNeuronposGaussian(this, varargin)
            % Initializes neuron position by sampling from a Gaussian
            % normal distribution standard deviation sigma and mean mu. 
            % By default, the mean value is in the middle of the muscle 
            % fibres and the standard deviation equals 0.3cm.
            this.System.f.initNeuronposGaussian(varargin{:});
        end
        
        function varargout = plot(this, varargin)
            % plot the surface EMG signal
            %
            % See also: EMGSystem
            [varargout{1:nargout}] = this.System.plot(varargin{:});
        end
        
        function h = plot1D(this, t, y, point)
            % plots the EMG signal at point over time
            geo = this.System.musclegeometry;
            if numel(point) ~= 3
                error('Input argument point must have 3 entries.\n');
            end
            if any([point(1)>geo(1), point(2)>geo(2), point(3)>geo(3)+geo(4)]) || any(point < 0)
                error('Input argument point must lie inside muscle.\n');
            end
            
            idx(1) = round(point(1)/this.System.h(1)) + 1; % ./ fails if point is row and h is column
            idx(2) = round(point(2)/this.System.h(2)) + 1;
            idx(3) = round(point(3)/this.System.h(3)) + 1;
            dim = this.System.dim;
            idx = idx(1) + dim(1)*(idx(2)-1) + dim(1)*dim(2)*(idx(3)-1);
            
            h = figure;
            plot(t, y(idx,:), t, abs(y(idx,:)),'LineWidth',2);
            title(sprintf('EMG signal at point (%g,%g,%g)',point(1),point(2),point(3)));
            xlabel('time [ms]');
            yl = '\phi_e';
            if point(3) > geo(3)
                yl = '\phi_o';
            end
            ylabel([yl, ' [mV]']);
            legend([yl, '(t)'], ['|',yl,'(t)|']);
        end
        
        function [reduced,time] = buildReducedModel(this, target_dim, DEIM_order)
            tic;
            if nargin < 3
                DEIM_order = find(this.Approx.SingularValues > 1e-7, 1, 'last');
                if nargin < 2
                    %target_dim = size(this.Data.V,2);
                    target_dim = find(this.SpaceReducer.SingularValues > 1e-4, 1, 'last');
                end
            end
            
            this.Approx.Order = DEIM_order;
            reduced = buildReducedModel@models.BaseFullModel(this, target_dim);
            time = toc;
        end
        
        function W = computeWfromV(this, V)%#ok
            % Idea 1: (Timm) use W = A*V
%             W = this.System.A.evaluate(1,0)*V;
            
            % Idea 2: Simply use zeros for the last row - it's a side
            % condition which is satisfied for the reduced basis anyways
            % (by definition as the values are generated from trajectories
            % that satisfy zero mean at each point, hence does the POD)
            W = [V; zeros(1,size(V,2))];
            
            o = general.Orthonormalizer;
            W = o.orthonormalize(W);
        end
        
        function off4_genApproximationTrainData(this)
            % Use W=V for off4_genApproximationTrainData, otherwise there 
            % is an error in data.ApproxTrainData.computeFrom
            % (because V and W don't have the same number of rows).
            oldW = this.Data.W;
            this.Data.W = this.Data.V;
            off4_genApproximationTrainData@models.BaseFullModel(this);
            this.Data.W = oldW;
        end
        
        function off5_computeApproximation(this)
            A = KerMor.App;
            host = A.getHost;
            
            
            % On LEAD: Write ApproxTrainData.fxi in a matrix - this makes
            % POD of fxi much faster
            if length(host) >= 4 && strcmp(host(1:4), 'lead')
                start = tic;
                atd = this.Data.ApproxTrainData;
                fxi = [];
                if ~isempty(atd.fxi)
                    for idx = 1:atd.fxi.getNumBlocks
                        fxi = [fxi atd.fxi.getBlock(idx)];
                    end
                end
                this.Approx.computeDEIM(this.System.f, fxi);
                this.OfflinePhaseTimes(5) = toc(start);
            else
                off5_computeApproximation@models.BaseFullModel(this);
            end
        end
    end
    
    methods(Static, Access=protected)
        function this = loadobj(this)
            if ~isa(this, 'EMGModel')
                sobj = this;
                this = EMGModel;
                this.RandSeed = sobj.RandSeed;
                this = loadobj@models.BaseFullModel(this, sobj);
            else
                this = loadobj@models.BaseFullModel(this);
            end
        end
    end
    
    methods  % plotting methods illustrating the motor unit configuration
        function plotMuscleConfiguration(this)
            this.plotMUcrosssection;
            mg = this.System.musclegeometry;
            fprintf('The muscle is %gcm long, %gcm wide and %gcm high.\n',mg(1),mg(2),mg(3));
            if (mg(4) ~= 0)
                fprintf('There is a %gcm high skin/ fat layer on top of the msucle.\n',mg(4));
            end
            h = this.System.h;
            fprintf('The discretization width is h = [%g, %g, %g]cm in x-, y- and z-direction.\n',h(1),h(2),h(3));
            fprintf('The muscle fibres are grouped into %g motor units.\n',this.System.f.numMU);
            if ~strcmp(this.System.f.FiringTimesMode,'precomputed')
                disp('The motoneurons'' firing times were precomputed by a motoneuron model.');
            elseif ~strcmp(this.System.f.FiringTimesMode,'simulate')
                disp('The motoneurons'' firing times are determined by simulation of a motoneuron model.');
            elseif ~strcmp(this.System.f.FiringTimesMode,'custom')
                disp('Custom motoneuron firing times are used.');
            end
            
        end
        
        function plotMU3d(this, numbers)
            % plots all fibres of one motor unit over the muscle.
            % number: index of the motor unit
            %rhs = this.System.f;
            sys = this.System;
            numstr = '';
            figure
            hold on
            if nargin < 2
                numbers = 1:min(this.System.f.numMU, 7);
            end
            ls = LineSpecIterator(length(numbers));
            for num = numbers
                pos=[];
                for i=1:sys.dim(2)
                    for j=1:sys.zdim_m
                        if this.MUdistribution(i,j) == num
                            pos = [pos,[i;j]];%#ok
                        end
                    end
                end
                for i = 1:size(pos,2)
                    X = [0,sys.musclegeometry(1)];
                    Y = [(pos(1,i)-1)*sys.h(2),(pos(1,i)-1)*sys.h(2)];
                    Z = [(pos(2,i)-1)*sys.h(3), (pos(2,i)-1)*sys.h(3)];
                    plot3(X, Y, Z, 'Color',ls.nextColor,'LineWidth',2)
                    plot3((this.neuronpos(pos(1,i),pos(2,i))-1)*sys.h(1),Y(1),Z(1),'xk','LineWidth',2)
                end
                numstr = [numstr, sprintf(', #%g',num)];%#ok
            end
            grid on
            
            title(['Muscle Fibres and Neuromuscular Junctions of Motor Unit(s) ',numstr]);
            xlabel('x-direction [cm]');
            ylabel('y-direction [cm]');
            zlabel('z-direction [cm]');
            axis([0 sys.musclegeometry(1) 0 sys.musclegeometry(2) ...
                0 sys.musclegeometry(3)+sys.musclegeometry(4)])
            view(-44,12);
            daspect([1,1,1]);
            hold off
        end
        
        function plotMUcrosssection(this, numbers)
            % plots the motor unit distribution over muscle cross section.
            % Musclefibres of the same motor unit have the same color and
            % marker.
            % Parameters:
            % circles: struct with fields centers (matrix 2 \times numMU) and radii (vector length numMU).
            
            rhs = this.System.f;
            if nargin < 2  % plot all motor units
                numbers = 1:rhs.numMU;
            end
            
            sys = this.System;
            geo = sys.musclegeometry;
            pos=cell(1,rhs.numMU);
            for i=1:sys.dim(2)
                for j=1:sys.zdim_m
                    pos{rhs.MUdistribution(i,j)}=[pos{rhs.MUdistribution(i,j)},[i;j]];
                end
            end
            figure
            % frames for body region
            if geo ~= 0
                X = [0, 0, geo(2), geo(2), 0];
                Y = [geo(3), geo(3)+geo(4), ...
                    geo(3)+ geo(4), geo(3), geo(3)];
                patch(X', Y', [0;0;0;0;0], 'FaceColor',[0.9 0.9 .9])
            end
            hold on
            
            ls = LineSpecIterator(length(numbers),1000);
            for num = numbers
                if isempty(pos{num})
                    continue;
                end
                col = ls.nextColor;
                plot((pos{num}(1,:)-1)*sys.h(2),(pos{num}(2,:)-1)*sys.h(3), '.', 'Marker', ls.nextMarkerStyle, 'MarkerEdgeColor', col, 'LineWidth', 1);
                if ~isempty(rhs.MUs{num}.centre) && rhs.numMU > 1
                    plotCircle(rhs.MUs{num}.centre,rhs.MUs{num}.radius);
                end
            end
            title('Motor Unit Distribution over Muscle Cross Section');
            xlabel('y-direction [cm]');
            ylabel('z-direction [cm]');
            axis([0,geo(2), 0,geo(3)+geo(4)])
            daspect([1 1 1]);
            hold off
            
            function plotCircle(cent, rad)
                dphi = 0.01;
                phi = 0:dphi:2*pi;
                x = cent(1) + rad*cos(phi);
                y = cent(2) + rad*sin(phi);
                if geo(4) > 0
                    idx = logical(y <= geo(3)*1.01);
                    x = x(idx);
                    y = y(idx);
                end
                plot(x,y,'-','Color',col);
            end
        end
    end
    
    methods  % getters and setters: interface to System and RHS
        function set.MUdistribution(this, value)
            if ~all(size(value) == [this.System.dim(2) this.System.zdim_m])
                error('Wrong dimension.')
            end
            
            numnew = max(max(value));
            % initialize MUs only if number of motor units changes.
            % otherwise, this corresponds to just rearranging the old MUs.
            if numnew ~= this.numMU  
                MUs = cell(1, numnew);
                types = rand(numnew); % RHS.rs is private:  types = this.System.f.rs.rand(numnew);
                for idx = 1:numnew
                    % ToDo check default values for MU
                    MUs{idx} = struct('type', types(idx), 'firing_times', [0], 'centre', [], 'radius', []);
                end
                this.System.f.MUs = MUs;
            end
            this.System.f.MUdistribution = value;
        end
        
        function val = get.MUdistribution(this)
            val = this.System.f.MUdistribution;
        end
        
        function value = get.numMU(this)
            value = max(max(this.MUdistribution));
        end
        
        function value = get.FiringTimesMode(this)
            value = this.System.f.FiringTimesMode;
        end
        
        function set.FiringTimesMode(this, value)
            this.System.f.FiringTimesMode = value;
        end
        
        function value = get.neuronpos(this)
            value = this.System.f.neuronpos;
        end
        
        function set.neuronpos(this, value)
            if ~all(size(value) == [this.System.dim(2) this.System.zdim_m])
                error('Wrong dimension.')
            end
            this.System.f.neuronpos = value;
        end

        function value = get.MUs(this)
            value = this.System.f.MUs;
        end
        
        function set.MUs(this, value)
            %if 
            %    error('')
            %end
            this.System.f.MUs = value;
        end
        
        function value = get.sigmao(this)
            value = this.System.sigmao;
        end
        
        function set.sigmao(this, value)
            if ~all(size(value) == [3 3]) || any(any(value-value')) 
               error('sigmao must be a diagonal 3 x 3 matrix.')
            end
            this.System.sigmao = value;
        end
        
        function value = get.sigmae(this)
            value = this.System.sigmae;
        end
        
        function set.sigmae(this, value)
            if ~all(size(value) == [3 3]) || any(any(value-value')) 
               error('sigmae must be a diagonal 3 x 3 matrix.')
            end
            this.System.sigmae = value;
        end
        
        function value = get.sigmai(this)
            value = this.System.sigmai;
        end
        
        function set.sigmai(this, value)
            if ~all(size(value) == [3 3]) || any(any(value-value')) 
               error('sigmai must be a diagonal 3 x 3 matrix.')
            end
            this.System.sigmai = value;
        end
        
%         function set.T(this, value)
%             this.T = T;
%             this.System.f.updateT;
%         end
    end
end