classdef Pool < handle
% Pool: A motorunit pool with either Fuglevand or Shorten model for force
% generation.
%
% @author Daniel Wirtz @date 2014-05-23
%
% @new{0,7,dw,2014-05-23} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(Dependent)
        % The fibre types.
        %
        % @default [] @type rowvec<double>
        % See also: FibreTypeWeights
        FibreTypes;
        
        % The fibre type weights.
        % Must be the same length as FibreTypes
        %
        % @default [] @type rowvec<double>
        % See also: FibreTypes
        FibreTypeWeights = [];
        
        % The number of fibres in the pool.
        %
        % @type integer @default 2
        NumFibres;
    end

    properties 
        % The type of the pool.
        %
        % Admissible values:
        % 1: Use the Fuglevand force model
        % 2: Use the Shorten force model
        %
        % 
        PoolType;
        
        % Flag for force normalization
        %
        % If true, the computed force will be normalized to 1.
        %
        % @default false @type boolean
        NormalizedForce = false;
    end
    
    properties(SetAccess=private)
        fugle;
        shorten;
        alpha;
        alpha_all;
    end
    
    properties(Access=private)
        lastdt;
        ft_ = [0 1];
        ftw_ = [.5 .5];
    end
    
    methods
        
        function this = Pool(type, shorten_model)
            % Creates a new motorunit pool
            %
            % Parameters:
            % type: The pool type. 1=Fuglevand, 2=Shorten @type integer
            % @default 1
            % shorten_model: The shorten model to use; only applies for
            % type==2. @type models.motorunit.Shorten @default Neumann
            % version
            if nargin < 2
                if nargin < 1
                    type = 1;
                end
                % We use the Neumann shorten version in the motorunit pools
                shorten_model = models.motorunit.Shorten('SarcoVersion',2);
                shorten_model.EnableTrajectoryCaching = true;
            end
            this.PoolType = type;
            this.fugle = models.motorunit.Fuglevand;
            this.shorten = shorten_model;
        end
        
        function [force, forces] = simulate(this, meancurrent, T, dt)
            % Simulates the pool for given mean current and times.
            %
            % Parameters:
            % meancurrent: The input mean current for motoneuron
            % activation. @type double
            % T: The final time `T` @type double
            % dt: The time step `\delta t` @type double
            %
            % Return values:
            % force: The weighted force over the times 0:dt:T. @type
            % rowvec<double>
            % forces: The single forces for each motorunit in the pool,
            % contained in each row. @type matrix<double>
            ft = this.ft_;
            nf = length(ft);
            t = 0:dt:T;
            forces = zeros(nf,length(t));
            this.fugle.Seed = 1;
            switch this.PoolType
                case 1
                    for fidx = 1:nf
                        forces(fidx,:) = this.fugle.getForces([ft(fidx); meancurrent],T,dt);
                    end
                case 2
                    m = this.shorten;
                    m.T = T;
                    m.dt = dt;
                    for fidx = 1:nf
                        
                        %MJ
                        % von 1 auf 2 gesetzt
                        [~,y] = this.shorten.simulate([ft(fidx); meancurrent],2);
                        %/MJ
                        
                        %Remove negative forces due to initial wobblyness
                        y(2,y(2,:)<0) = 0;
                        forces(fidx,:) = y(2,:);
                    end
            end
            
            % Weighted sum of all forces is the total force
            force = this.ftw_*forces;
            
            if this.NormalizedForce
                mval = max(force);
                if mval > 0
                    forces = forces/mval;
                    force = force/mval;
                end
            end
        end
        
        function prepare(this, meancurrent, T, dt)
            % Simulates the pool for given mean current and times and
            % stores the result for later access via
            % models.motorunit.Pool#getActivation.
            %
            % Parameters:
            % meancurrent: The input mean current for motoneuron
            % activation. @type double
            % T: The final time `T` @type double
            % dt: The time step `\delta t` @type double
            %
            % See also: getActivation
            if KerMor.App.Verbose > 0
                fprintf('Preparing motorunit pool for mean current %g over times 0:%g:%g\n',meancurrent,dt,T);
            end
            % Increase the internal resolution - live muscle simulations
            % with same accuracy as desired output often require finer
            % internal time steps.
            dt = max(.001, dt / 100);
            [this.alpha, this.alpha_all] = this.simulate(meancurrent, T, dt);
            this.lastdt = dt;
        end
        
        function [alpha, alpha_all]  = getActivation(this, t)
            % Returns the activation level for a given time instance after
            % use of models.motorunit.Pool#prepare
            %
            % Return values:
            % force: The weighted force over the times 0:dt:T. @type
            % rowvec<double>
            %
            % See also: prepare
            alpha = [];
            alpha_all = [];
            a = this.alpha;
            if ~isempty(a)
                pos = min(size(a,2),1+round(t/this.lastdt));
                alpha = a(pos);
                alpha_all = this.alpha_all(:,pos);
            end
        end
        
        function value = get.FibreTypes(this)
            value = this.ft_;
        end
        
        function value = get.FibreTypeWeights(this)
            value = this.ftw_;
        end
        
        function set.FibreTypeWeights(this, value)
            % Ensure its a row vector
            value = reshape(value,1,[]);
            if abs(sum(value)-1) > 1e-8
                error('FibreTypeWeights must add up to one!');
            end
            if ~isempty(this.ft_) && length(this.ft_) ~= length(value)
                error('FibreTypeWeights must have the same length as FibreTypes!')
            end
            this.ftw_ = value;
        end
        
        function set.FibreTypes(this, value)
            % Ensure its a row vector
            value = reshape(value,1,[]);
            % Ensure we have a compatible length
            if length(this.ftw_) ~= length(value)
                this.ftw_ = ones(size(value))/length(value);
            end
            this.ft_ = value;
        end
    end
    
    methods(Static)
        
        function res = test_MotorunitPool
            res = true;
            
            p = models.motorunit.Pool;
            p.prepare(4,10,.1);
            
            p = models.motorunit.Pool(1);
            p.FibreTypes = rand(1,50);
            v = rand(1,50);
            p.FibreTypeWeights = v / sum(v);
            p.NormalizedForce = true;
            [f,allf] = p.simulate(4,500,.1);
            figure;
            t = 0:.1:500;
            plot(t,f,'r',t,allf,'g');
            % Need max force equal to one!
            res = res && abs(max(f)-1)<1e-8;
            
            p = models.motorunit.Pool(2);
            p.FibreTypes = rand(1,4);
            v = rand(1,4);
            p.FibreTypeWeights = v / sum(v);
            p.prepare(9,20,.01);
            
            m = models.motorunit.Shorten('SarcoVersion',1);
            p = models.motorunit.Pool(2,m);
            T = 60;
            [f,allf] = p.simulate(4,T,.01);
            p.prepare(4,T,.01);
            t = 0:.01:T;
            a = p.getActivation(t);
            res = res && max(f(:)-a(:)) < 0.002;
            figure;
            plot(t,f,'r',t,allf,'g');
        end
        
    end
    
end