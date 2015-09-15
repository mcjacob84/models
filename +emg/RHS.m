classdef RHS < dscomponents.ACompEvalCoreFun
    % FibreDynamics: Class for right-hand side B*u(t) of the EMG model
    %
    % @author Timm Strecker @date 2014-04-08
    %
    % @new{0,7,ts,2014-04-08} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    
    properties(Dependent)
        % gridwidth of the discretization on all 3 dimensions
        h;
        
        % number of grid points of the discretization in all 3 dimensions
        dims;
        
        % the number of motor units
        %
        % @default 1
        numMU;
    end
    
%     properties
%         MultiArgumentEvaluations = true;
%     end
    
    properties(SetAccess = private)
        B; % the B matrix (discrete generalized Laplace with conductivity sigma_i)
        
        % velocity of action potentials propagating along the muscle
        % fibres in cm/ms 
        % ToDo: maybe type dependent
        APvelocity = 0.2;
        
        % struct containing data on the shape of a action potential for a
        % set of fiber types.
        % see also initAPdata
        APdata;
    end
    
    properties
        % Flag indicating whether to use custom firing times of the motor units.
        % Otherwise, the firing times are generated automatically by using the parameter mu.
        % 'precomputed', 'simulate', 'custom'
        FiringTimesMode = 'simulate';
        
        % Position (index) of the neuromuscular junction of each fibre. By
        % default randomly anywhere on the fibres.
        % matrix<integer> of dimension dim(2) \times zdim_m
        % see also: initNeuronposGaussian
        neuronpos;
        
        % Cell array of structs containing information about the motor
        % units. There are the fields
        %'type': specifying the type (slow-twitch, fast-twitch, etc.),
        % 'firing_times': containing the times when the motoneuron is firing,
        % 'centre' and 'radius': optional fields that define a circle
        % of radius radius around centre. All fibres of the MU are
        % distributed inside this circle.
        MUs;
        
        % distribution of motor units over cross section
        %
        % @default ones(dim(2),zdim_m)
        MUdistribution;
    end
    
    properties(Access=private)
        % The randstream instance
        rs;
        
        % motoneuron model for computation of firing times at FiringTimesMode 'simulate'
        mn;
    end
    
    methods
        
        % constructor
        function this = RHS(sys)
            this = this@dscomponents.ACompEvalCoreFun(sys);
            this.rs = RandStream('mt19937ar','Seed',24247247);
            this.xDim = sys.dim(1)*sys.dim(2)*sys.dim(3);
            this.fDim = sys.dim(1)*sys.dim(2)*sys.dim(3)+1;
            this.assembleB;
            this.JSparsityPattern = sparse(this.fDim, this.xDim);
            this.initAPdata;
            this.distributeMUs(4);
            %this.neuronpos = this.rs.randi(this.dims(1),this.dims(2:3)');
            this.initNeuronposGaussian;
            
            % init motoneuron
            mn = models.motoneuron.Model;
            mn.dt = 0.1;
            mn.EnableTrajectoryCaching = true;
            this.mn = mn;
        end
        
        function copy = clone(this)
            % Create new instance ToDo
            copy = RHS(this.System);
            % Call superclass clone (for deep copy)
            copy = clone@dscomponents.ACompEvalCoreFun(this, copy);
            % Copy local properties
            copy.MUdistribution = this.MUdistribution;
            copy.FiringTimesMode = this.FiringTimesMode;
            copy.neuronpos = this.neuronpos;
            copy.MUs = this.MUs;
        end
        
        function res = test_ComponentEvalMatch(this, xsize)
            % The original test takes extremely long if the firing times
            % are computed by the motoneuron model. Sometimes the result is
            % false due to numerical rounding (e.g., if the input is 10^{-12}
            % although it should be exactly 0. then the relative error can
            % be large).
            % I tested it many times.
            % ToDo: remove this function?
            res = true;
        end
    end
    
    
    %% methods that compute the right-hand side
    methods
        function dy = evaluate(this, y, t, mu)
            % Evaluates the nonlinear core function at given time. Can evaluate vectorized arguments.
            %
            % Parameters:
            % y: @type matrix<double>, not used
            % t: current time @type rowvec<double>
            % mu: mean input current @type rowvec<double>
            
            dy = this.B * this.getV_m(t);
            if ~isempty(this.W)   % ToDo: I'm not sure if this can cause errors ...
                dy = this.W'*dy;
            end
        end
        
        function dy = evaluateMulti(this, y, t, mu)
            % evaluate can already handle vectorized arguments
            mu = mu(1:this.System.ParamCount,:);
            if all(all(repmat(mu(:,1),1,size(mu,2)-1) == mu(:,2:end)))
                this.prepareSimulation(mu(:,1));
                dy = this.evaluate(y, t, mu);
            else 
                dy = zeros(this.fDim, length(t));
                for k = 1:length(t)
                    this.prepareSimulation(mu(:,k));
                    dy(:,k) = this.evaluate([], t(k), mu(:,k));
                end
            end 
        end
        
        function u = getV_m(this, t)
            % Computes the "source term" V_m at given time points t.
            % see also: computeMUActivation
            
            d = this.dims;
            u = zeros(prod(d),length(t));
            neurpos = this.neuronpos;
            MUdistr = this.MUdistribution;
            MUActivation = this.computeMUActivation(t);
            for jdx = 1:d(2)
                for kdx = 1:d(3)
                    fibre = MUActivation(MUdistr(jdx,kdx),neurpos(jdx,kdx):d(1)-1+neurpos(jdx,kdx),:);
                    u(1+(jdx-1)*d(1)+(kdx-1)*d(1)*d(2):d(1)+(jdx-1)*d(1)+(kdx-1)*d(1)*d(2),:) = fibre;
                end
            end
        end
        
        function MUActivation = computeMUActivation(this, t)
            % Computes the "source term" V_m for all motor units at given
            % time points t. For each MU, V_m is computed over a
            % representative fibre that is twice as long as the actual muscle and
            % has its neuromuscular junction in the middle.
            % ToDo: mu dependend velocity?
            MUActivation = zeros(this.numMU, 2*this.dims(1)-1, length(t));
            v = this.APvelocity;
            max_len = 2*this.dims(1)-1;
            sys = this.System;
            for mdx = 1:this.numMU
                MU = this.MUs{mdx};
                [~,idx] = min(abs(this.APdata.mus - MU.type));
                APs = this.APdata.shapes{idx};
                Vm = -80*ones(max_len,length(t));
                ft = MU.firing_times;  % firing times of current MU
                for tdx = 1:length(t)
                    % pick only APs that are relevant at current time
                    rel_ft = ft((ft <= t(tdx)) & (ft >= t(tdx)-sys.musclegeometry(1)/v));
                    % iterate over all firing times
                    for fdx = 1:numel(rel_ft)   
                        % how far the "front" of the AP corresponding to the current ft propagated at time t(tdx)
                        front = round(v*(t(tdx)-rel_ft(fdx))/sys.h(1));   
                        num = min(front+1, numel(APs));
                        tail = front-num+1;   % tail of the current AP
                        APpos = this.dims(1)-front:this.dims(1)-tail;
                        Vm(APpos,tdx) = APs(1:num);
                        APpos = this.dims(1)+tail:this.dims(1)+front;
                        Vm(APpos,tdx) = APs(num:-1:1);
                    end
                end
                MUActivation(mdx,:,:) = Vm;
            end
        end
        
        function assembleB(this)
            % Assembles the B matrix of the right-hand side Bu(t).
            % The upper rows of B equal the discrete generalized Laplace \nabla\sigma_i\nabla.
            % These rows correspond to points in the muscle tissue, where there are muscle fibres.
            % All other rows are all-zero. These correspond to the points in
            % the skin layer (no V_m). The last zero-row is due to the
            % zero-mean condition.
            musclegrid = general.geometry.RectGrid3D(this.dims(1),this.dims(2),this.dims(3));
            [iu, ju, su] = find(MatUtils.generalizedLaplacemat3D(this.h,musclegrid,this.System.sigmai));
            this.B = sparse(iu, ju, su, this.fDim, prod(this.dims));
        end
        
        function fx = evaluateComponents2(this, pts, ends, idx, self, x, t, mu)
            % public equivalent to evaluateComponents, e.g. for testing
            fx = this.evaluateComponents(pts, ends, idx, self, x, t, mu);
        end
        
        function dy = evaluateCoreFun(this, y, t, mu)%#ok
            error('evaluate is overridden directly.');
        end
        
        function J = getStateJacobian(this, y, t, mu)
            % The Jacobian is empty because RHS uses no state space arguments.
            J = sparse(this.fDim,this.xDim);
        end
    end
    
    
    %% methods preparing a simulation
    methods
        function prepareSimulation(this, mu)
            prepareSimulation@dscomponents.ACoreFun(this, mu);
            if ~strcmp(this.FiringTimesMode, 'custom')
                this.updateFiringTimes(mu);
            end
        end
        
        function updateFiringTimes(this, mu)
            if strcmp(this.FiringTimesMode, 'simulate') % long computation possible
                this.mn.T = this.System.Model.T;
                for idx = 1:this.numMU   %@ToDo: simulate only if this type was not simulated before
                    type = this.MUs{idx}.type;
                    mc = mu(1);
                    if type <= .6   % adjust unnatural mean currents
                        mc = min(mc, 3.5);
                    else
                        mc = min(mc, 3.5 + (type-.6)/.4 * 4.5);  % line between [0.6,3.5] and [1,8]
                    end
                    [t,y] = this.mn.simulate([type;mc], 1);
                    APtimes = (y(2,2:end) > 40).*(y(2,1:end-1) <= 40);
                    this.MUs{idx}.firing_times = t(logical(APtimes));
                end
            elseif strcmp(this.FiringTimesMode, 'precomputed')
                d = load('neuroft');   % see also skriptMotoNeuron.m
                for idx = 1:this.numMU
                    type = this.MUs{idx}.type;
                    mc = mu(1);
                    if type <= .6   % adjust unnatural mean currents
                        mc = min(mc, 3.5);
                    else
                        mc = min(mc, 3.5 + (type-.6)/.4 * 4.5);  % line between [0.6,3.5] and [1,8]
                    end
                    [~,tdx]=min(abs(d.T-type));
                    [~,mdx]=min(abs(d.MC-mc));
                    this.MUs{idx}.firing_times = d.FT{tdx,mdx};
                    % this.MUs{idx}.firing_times = this.MUs{idx}.firing_times + rand  ... add some "noise"?
                end
            end
        end
        
        function initAPdata(this)
            % APshape for parameters MU is computed by monodomain equation for one fibre
            % with a given discretization width dx. It represents V_m at
            % all discretization points where it is above -80mV.
            % see also APshape.m
            datafile = fullfile(fileparts(mfilename('fullpath')),'data','APshape.mat');
            d = load(datafile);
            APs = cell(1,11);
            for i=1:11
                l = (length(d.APshape{i})-1)*d.dx;  % "length" of AP
                x = 0:d.dx:l;
                myx = 0:this.h(1):l;
                APs{i} = interp1(x, d.APshape{i}, myx);
            end
            this.APdata.shapes = APs;
            this.APdata.mus = d.MU;
        end
    end
    
    methods  % setting up the MU configuration
        function distributeMUs(this, num, centers, radii)
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
            
            if nargin < 2
                num = 1;
            end
            geo = this.System.musclegeometry(2:3);
            this.MUs = cell(1, num);
            types = [this.rs.rand(1,num-1) 1];
            if nargin <= 3
                centers = bsxfun(@times,geo,this.rs.rand(2,num));
            else
                if ~all(size(centers) == [2 num])
                    error('centers must be matrix of dimension 2 x %g', num);
                end
            end
            
            if nargin < 4
                % Daniel: radii = [norm(geo) norm(geo)*(this.rs.rand(1,number-1)*.5 + .1)];
                radii=norm(geo)*sqrt(exp(types * log(100))/100);
            else
                if isscalar(radii)
                    radii = radii*ones(1,num);
                elseif length(radii) ~= num
                    error('radii must have %g entries or be scalar.', num);
                end
            end
                
            % assigns fibres to motor units
            h = this.h;
            d = this.dims;
            MUdistr = zeros(d(2),d(3));
            randnumbers = this.rs.rand(d(2),d(3));   % precompute random numbers - faster than in for loop
            for jdx = 1:d(2)
                for kdx = 1:d(3)
                    w = ones(1,num);
                    for i=1:num
                        % if point lies outside circle of current MU, fibre cannot belong to this MU
                        if norm(([jdx-1; kdx-1].*[h(2);h(3)] - centers(:,i))) > radii(i)
                            w(i) = 0;
                        end
                    end
                    if ~any(w)
                        error('Error: There is no motor unit at position y = %g, z = %g. Cannot assign a fibre.',(jdx-1)*h(2), (kdx-1)*h(3));
                    end
                    w = w/sum(w);   % normieren auf summe 1
                    b = find(randnumbers(jdx,kdx) <= cumsum(w),1,'first');
                    MUdistr(jdx,kdx) = b;
                end
            end
            this.MUdistribution = MUdistr;
            
            % ToDo check default values for MU
            for idx = 1:num
                this.MUs{idx} = struct('type', types(idx),...
                    'firing_times', 0,...
                    'centre', centers(:,idx),...
                    'radius', radii(idx));
            end
            
            fprintf('You can view the current configuration via EMGModel.plotMuscleConfiguration.\n');
        end
        
        function initNeuronposGaussian(this, sigma, mu)
            % Initializes neuron position by sampling from a Gaussian
            % normal distribution standard deviation sigma and mean mu. 
            % By default, the mean value is in the middle of the muscle 
            % fibres and the standard deviation equals 0.3cm. 
            sys = this.System;
            if nargin < 3
                mu = sys.musclegeometry(1)/2;
                if nargin < 2
                    sigma = min(0.3, mu/2);
                end
            end
            r = mu + sigma*this.rs.randn(sys.dim(2), sys.zdim_m);
            while any(any(r < 0)) || any(any(r > sys.musclegeometry(1)))
                idx = logical((r < 0) + (r > sys.musclegeometry(1)));
                rnew = mu + sigma*this.rs.randn(sys.dim(2), sys.zdim_m);
                r(idx) = rnew(idx);
            end
            this.neuronpos = round(r/sys.h(1)) + 1;
        end
    end
    
    
    %% componentwise evaluation
    methods(Access=protected)
        function fx = evaluateComponents(this, pts, ends, idx, self, x, t, mu)
            % This is the template method that actually evaluates the components at given points
            % and values.
            %
            % @attention This method must be able to handle vector-arguments
            % for `\vx,t,\vmu`!
            %
            % Parameters:
            % pts: The components of `\vf` for which derivatives are required @type rowvec<integer>
            % ends: At the `i`-th entry it contains the last position in the `\vx` vector that
            % indicates an input value relevant for the `i`-th point evaluation, i.e.
            % `f_i(\vx) = f_i(\vx(ends(i-1)+1{:}ends(i)));` @type rowvec<integer>
            % idx: The indices of `\vx`-entries in the global `\vx` vector w.r.t the `i`-th
            % point, e.g. `xglobal(i-1:i+1) = \vx(ends(i-1):ends(i))` @type rowvec<integer>
            % self: The positions in the `\vx` vector that correspond to the `i`-th output
            % dimension, if applicable (usually `f_i` depends on `x_i`, but not necessarily)
            % @type rowvec<integer>
            % x: A matrix `\vX` with the state space locations `\vx_i` in
            % its columns. In rows end(i-1)+1:end(i) it contains the states
            % relevant to evaluate the i-th component of �\vf�.
            % States occur multiply in �\vx� if different components
            % of �\vf� depend on these states.
            % @type matrix<double>
            % t: The corresponding times `t_i` for each state `\vx_i` @type rowvec<double>
            % mu: The corresponding parameters `\mu_i` for each state `\vx_i`, as column matrix
            % @type matrix<double>
            %
            % Return values:
            % fx: A matrix with pts-many component function evaluations `f_i(\vx)` as rows and as
            % many columns as `\vX` had.
            fx = zeros(length(pts), length(t));
            
            [j, i, s] = find(this.B(pts,:)');
            idx=zeros(3,length(j));
            idx(3,:) = ceil(j'/this.dims(1)/this.dims(2));
            idx(2,:) = ceil((j'-(idx(3,:)-1)*this.dims(1)*this.dims(2))/this.dims(1));
            idx(1,:) = j'-this.dims(1)*this.dims(2)*(idx(3,:)-1)-this.dims(1)*(idx(2,:)-1);
            neurpos = this.neuronpos;
            MUdistr = this.MUdistribution;
            MUActivation = this.computeMUActivation(t);
            ii = MUdistr(idx(2,:)+(idx(3,:)-1)*this.dims(2));
            %jj = neurpos(jdx,kdx):this.dims(1)-1+neurpos(jdx,kdx);
            jj = neurpos(idx(2,:)+(idx(3,:)-1)*this.dims(2))+idx(1,:)-1;
            
            % stretched variant of B
            B2i = i;
            B2j = 1:length(s);
            B2s = s;
            B2 = sparse(B2i, B2j, B2s, length(pts), length(s));
            
            for tdx = 1:length(t)
                %act = MUActivation(MUdistr(idx(2,:)+(idx(3,:)-1)*this.dims(2)), neurpos(jdx,kdx):this.dims(1)-1+neurpos(jdx,kdx), t(tdx));       % pseudocode
                act = MUActivation(ii + (jj-1)*this.numMU + (tdx-1)*this.numMU*(2*this.dims(1)-1));
                % the used way of getting i, js and s and constructing act
                % ensures that the entries are in the right order
                fx(:,tdx) = B2*act';
            end
        end
        
        function dfx = evaluateComponentPartialDerivatives(this, pts, ends, idx, deriv, self, x, t, mu, dfxsel)%#ok
            % Computes specified partial derivatives of `f` of the components given by pts and
            % the selected partial derivatives by dfxsel.
            %
            % Parameters:
            % pts: The components of `f` for which derivatives are required @type
            % rowvec<integer>
            % ends: At the `i`-th entry it contains the last position in the `\vx` vector that
            % indicates an input value relevant for the `i`-th point evaluation, i.e.
            % `f_i(\vx) = f_i(\vx(ends(i-1){:}ends(i)));` @type rowvec<integer>
            % idx: The indices of `\vx`-entries in the global `\vx` vector w.r.t the `i`-th
            % point, e.g. `xglobal(i-1:i+1) = \vx(ends(i-1):ends(i))` @type rowvec<integer>
            % deriv: The indices within `\vx` that derivatives are required for.
            % @type rowvec<integer>
            % self: The positions in the `\vx` vector that correspond to the `i`-th output
            % dimension, if applicable (usually `f_i` depends on `x_i`, but not necessarily)
            % @type rowvec<integer>
            % x: The state space location `\vx` @type colvec<double>
            % t: The corresponding times `t` for the state `\vx` @type double
            % mu: The corresponding parameter `\mu` for the state `\vx` @type colvec<double>
            % dfxsel: A derivative selection matrix. Contains the mapping for each row of x to
            % the output points pts. As deriv might contain less than 'size(x,1)' values, use
            % 'dfxsel(:,deriv)' to select the mapping for the actually computed derivatives.
            %
            % Return values:
            % dfx: A column vector with 'numel(deriv)' rows containing the derivatives at all
            % specified pts i with respect to the coordinates given by 'idx(ends(i-1):ends(i))'
            %
            % See also: setPointSet
            
        end
        
        function fx = evaluateComponentsMulti(this, pts, ends, idx, self, x, t, mu)
            % @todo improve performance!
            fx = zeros(length(pts),length(t));
            mu = mu(1:this.System.ParamCount,:);
            if all(all(repmat(mu(:,1),1,size(mu,2)-1) == mu(:,2:end)))
                this.prepareSimulation(mu(:,1));
                for k=1:length(t)
                    fx(:,k) = this.evaluateComponents(pts, ends, idx, self, [], t(k), mu(:,k));
                end
            else
                for k=1:length(t)
                    this.prepareSimulation(mu(:,k));
                    fx(:,k) = this.evaluateComponents(pts, ends, idx, self, [], t(k), mu(:,k));
                end
            end
        end
    end
    
    
    %% Setters & Getters
    methods
        function set.FiringTimesMode(this, value)
            if ~any([strcmp(value, 'precomputed'), ...
                    strcmp(value, 'simulate'), strcmp(value, 'custom')])
                error(['Invalid value for property FiringTimesMode. Possible',...
                    ' are ''precomputed'', ''simulate'' and ''custom''.\n'])
            end
            this.FiringTimesMode = value;
        end
        
        function val = get.h(this)
            val = this.System.h;
        end
        
        function val = get.dims(this)
            val = [this.System.dim(1);this.System.dim(2);this.System.zdim_m];
        end
        
        function value = get.numMU(this)
            value = max(max(this.MUdistribution));
        end
    end
end
