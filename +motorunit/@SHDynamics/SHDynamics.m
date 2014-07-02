classdef SHDynamics < dscomponents.ACoreFun
    % SHDynamics: Class for nonlinear dynamics of shorten motorunit model
    %
    % @author Daniel Wirtz @date 2012-11-22
    %
    % @new{0,7,dw,2012-11-22} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    
    properties
        % The `V_s` value of the motoneuron input at which the MSLink_MaxFactor should be attained
        %
        % @type double @default 40
        MSLink_MaxFactorSignal = 40;
        
        % The maximal factor with which the `V_s` value of the motoneuron should be amplified
        % when added to the sarcomere equations
        %
        % @type double @default 7
        MSLink_MaxFactor = 7;
        
        % The minimal factor at which the `V_s` value of the motoneuron should be amplified
        % when added to the sarcomere equations.
        %
        % Between both limits, an exponential Gaussian weight is applied with a certain radius,
        % so that the amplification factor has its minimum value at zero and maximum value at
        % MSLink_MaxFactorSignal.
        %
        % See also FibreDynamics.getLinkFactor
        %
        % @type double @default .3
        MSLink_MinFactor = .3;
    end
    
    properties(SetAccess=private)
        % constants for sarcomere submodel specific for slow twitch fibre
        %
        % See also: Dynamics.initSarcoConst
        SarcoConst_slow;
        
        % constants for sarcomere submodel specific for fast twitch fibre
        %
        % See also: Dynamics.initSarcoConst
        SarcoConst_fast;
        
        % basis constant set for sarcomere submodel
        %
        % See also: Dynamics.initSarcoConst
        SarcoConst_base;
        
        % Positions of type-dependent constants
        %
        % See also: Dynamics.initSarcoConst
        SarcoConst_dynpos;
    end
    
    properties(Access=private)
        spm;
    end
    
    properties(Transient)
        hadpeak = false;
        motoconst;
        sarcoconst;
        signalweight;
    end
    
    methods
        function this = SHDynamics(sys)
            this = this@dscomponents.ACoreFun(sys);
            
            this.initSarcoConst;
            this.initJSparsityPattern;
            this.fDim = sys.dsa + sys.dm;
            this.xDim = this.fDim;
            this.spm = sys.Model.SinglePeakMode;
        end
        
        function prepareSimulation(this, mu)
            prepareSimulation@dscomponents.ACoreFun(this, mu);
            
            % Precompute motoneuron/sarcomere constants
            this.motoconst = this.getMotoConst(mu(1));
            this.sarcoconst = this.getSarcoConst(mu(1));
            
            this.hadpeak = false;
        end
        
        function dy = evaluate(this, y, t)
            % Evaluates the nonlinear core function at given time and
            % state. Can evaluate vectorized arguments. Here, each column
            % represents one state.
            % The parts for motoneuron, spindle and sarcomeres are evaluated seperately in
            % respective rate functions and the rates are concatenated in this function. The sarcomere part is reshaped such
            % that a matrix with 58 rows and one column for each single sarcomere is passed to the function
            % for (single) sarcomere dynamics.
            %
            % Parameters:
            % y: @type matrix<double>, each column is one state vector
            % t: current time @type rowvec<double>
            % mu: fibre type parameter @type rowvec<double> 0 = slow twitch fibre, 1 = fast twitch fibre
            
            dm = this.System.dm;
            
            %% Inner dynamics (with partial nonlinear linking)
            dy = zeros(size(y));
            
            %% Neuro dynamics (last 2 arguments are 2 phases: phase_prim, phase_sec)
            c = this.motoconst;
            % dendrites
            dy(1) = (-c(1)*(y(1)-c(11))-c(5)*(y(1)-y(2)))/c(7);
            % soma
            dy(2) = (-c(6)*(y(2)-c(11))-c(5)*(y(2)-y(1))...
                -c(4)*y(3).^3*y(4)*(y(2)-c(9))...
                -c(2)*y(5).^4*(y(2)-c(10))...
                -c(3)*y(6).^2*(y(2)-c(10)))/c(8);
            % the four gating variables
            dy(3) = 0.32*(13-y(2))/(exp((13-y(2))/5)-1)*(1-y(3))...
                -0.28*(y(2)-40)/(exp((y(2)-40)/5)-1)*(y(3));
            dy(4) = 0.128*(exp((17-y(2))/18))*(1-y(4))-4/(exp((40-y(2))/5)+1)*(y(4));
            dy(5) = 0.032*(15-y(2))/(exp((15-y(2))/5)-1)*(1-y(5))...
                -0.5*(exp((10-y(2))/40))*(y(5));
            dy(6) = 3.5/(exp((55-y(2))/4)+1)*(1-y(6))-0.025*(y(6));
            
            %% Sarcomere dynamics
            dy(dm+1:end) = this.dydt_sarcomere(y(dm+1:end),t);
            
            %% Link of motoneuron to sarcomer cell
            if this.spm && this.hadpeak
                signal = 0;
            else
                fac = this.MSLink_MaxFactor;
                if y(2) < this.MSLink_MaxFactorSignal
                    fac = this.getLinkFactor(y(2));
                end
                % "Link" at the only sarcomere
                signal = fac*y(2)/this.sarcoconst(1);
            end
            
            dy(dm+1,:) = dy(dm+1,:) + signal;
        end
        
        function dy = evaluateCoreFun(this, y, t, mu)%#ok
            error('evaluate is overridden directly.');
        end
        
        function J = getStateJacobian(this, y, t)
            % Evaluates the full Jacobian of the nonlinear core function.
            %
            % Parameters:
            % y: @type matrix<double>, each column is one state vector
            % t: current time @type rowvec<double>
            % mu: fibre type parameter, 0 = slow twitch fibre, 1 = fast twitch fibre, @type rowvec<double>
           
            s = this.System;
            dm = s.dm;
            dsa = s.dsa;
            bs = dm + dsa;
            J = spalloc(bs,bs,45+303);  % 45: 40 entries in motoneuron and spindle part + reserve, 303: entries per single sarcomere
            
            %% Neuron jacobian
            c = this.motoconst;
            JN = zeros(6,6);
            JN(1,1) = -(c(1) + c(5))/c(7);
            JN(2,1) = c(5)/c(8);
            JN(1,2) = c(5)/c(7);
            JN(2,2) = -(c(4)*y(4)*y(3)^3 + c(2)*y(5)^4 + c(3)*y(6)^2 + c(5) + c(6))/c(8);
            JN(3,2) = (8*(y(3) - 1))/(25*(exp(13/5 - y(2)/5) - 1)) - (7*y(3))/(25*(exp(y(2)/5 - 8) - 1)) + (y(3)*exp(y(2)/5 - 8)*((7*y(2))/25 - 56/5))/(5*(exp(y(2)/5 - 8) - 1)^2) + (exp(13/5 - y(2)/5)*((8*y(2))/25 - 104/25)*(y(3) - 1))/(5*(exp(13/5 - y(2)/5) - 1)^2);
            JN(4,2) = (8*exp(17/18 - y(2)/18)*(y(4) - 1))/1125 - (4*y(4)*exp(8 - y(2)/5))/(5*(exp(8 - y(2)/5) + 1)^2);
            JN(5,2) = (y(5)*exp(1/4 - y(2)/40))/80 + (4*(y(5) - 1))/(125*(exp(3 - y(2)/5) - 1)) + (exp(3 - y(2)/5)*((4*y(2))/125 - 12/25)*(y(5) - 1))/(5*(exp(3 - y(2)/5) - 1)^2);
            JN(6,2) = -(7*exp(55/4 - y(2)/4)*(y(6) - 1))/(8*(exp(55/4 - y(2)/4) + 1)^2);
            JN(2,3) = (3*c(4)*y(3)^2*y(4)*(c(9) - y(2)))/c(8);
            JN(3,3) =  ((8*y(2))/25 - 104/25)/(exp(13/5 - y(2)/5) - 1) - ((7*y(2))/25 - 56/5)/(exp(y(2)/5 - 8) - 1);
            JN(2,4) = (c(4)*y(3)^3*(c(9) - y(2)))/c(8);
            JN(4,4) =  - (16*exp(17/18 - y(2)/18))/125 - 4/(exp(8 - y(2)/5) + 1);
            JN(2,5) = (4*c(2)*y(5)^3*(c(10) - y(2)))/c(8);
            JN(5,5) = ((4*y(2))/125 - 12/25)/(exp(3 - y(2)/5) - 1) - exp(1/4 - y(2)/40)/2;
            JN(2,6) = (2*c(3)*y(6)*(c(10) - y(2)))/c(8);
            JN(6,6) = - 7/(2*(exp(55/4 - y(2)/4) + 1)) - 1/40;
            J(1:6,1:6) = JN;
            
            %% Sarcomere jacobian
            J(dm+1:bs,dm+1:bs) = this.Jac_Sarco(y(dm+1:end), t);
            fac = this.MSLink_MaxFactor;
            if (y(2) < this.MSLink_MaxFactorSignal)
                fac = this.getLinkFactor(y(2));
            end
            % link position is dm+1
            J(dm+1,2) = fac/this.sarcoconst(1);
        end
        
        function copy = clone(this)
            % Create new instance
            copy = models.motorunit.Dynamics(this.System);
            
            % Call superclass clone (for deep copy)
            copy = clone@dscomponents.ACompEvalCoreFun(this, copy);
            
            % Copy local properties
            % No local properties are to be copied here, as so far
            % everything is done in the constructor.
        end
        
        function initJSparsityPattern(this)
            % Initializes the Sparsity pattern of the Jacobian
            %
            s = this.System;
            dm = s.dm; dsa = s.dsa;
            size = dm + dsa;
            J = sparse(size,size);
            
            % Neuro
            i = [1,1,2,2,2,2,2,2,3,3,4,4,5,5,6,6];
            j = [1,2,1,2,3,4,5,6,2,3,2,4,2,5,2,6];
            J(1:dm,1:dm) = sparse(i,j,true,6,6);
            
            % Sarco
            mc = metaclass(this);
            s = load(fullfile(fileparts(which(mc.Name)),'JSparsityPattern'));
            J(dm+1:size,dm+1:size) = s.JP;
            
            linkpos = dm+1;
            J(linkpos,2) = true;
            this.JSparsityPattern = logical(J);
        end
        
        function f = getLinkFactor(this, y)
            % scalar link factor motoneuron to middle sarcomere
            %
            % See also: FibreDynamics.evaluate
            % Parameters:
            % y: @type rowvec<double>  second state from motoneuron
            % submodel
            f = this.MSLink_MinFactor + exp(-(y-this.MSLink_MaxFactorSignal).^2/150)...
                *(this.MSLink_MaxFactor-this.MSLink_MinFactor);
        end
    end
    
    methods(Access=private)
        function sc = getSarcoConst(this, mu_fibretype)
            % Precompute sacromere constants
            sc = this.SarcoConst_base;
            sc(this.SarcoConst_dynpos) = (1-mu_fibretype(1))*this.SarcoConst_slow + mu_fibretype(1)*this.SarcoConst_fast;
        end
        
        function c = getMotoConst(~, mu_fibretype)
            % getMotoConst: private getter function for motoneuron
            % constants, vectorial implementation
            %
            % mu_fibretype values are assumed to be in [0,1]
            
            % Membrane capacitance (NOT to be confused with Cm of the
            % sarcomere model!)
            Cm=1;
            
            Ri=70/1000;
            c = zeros(11,length(mu_fibretype));
            % cf. Cisi and Kohn 2008, Table 2, page 7
            Rmd = 14.4+6.05-coolExp(6.05,14.4,mu_fibretype);
            Rms=1.15+0.65-coolExp(0.65,1.15,mu_fibretype);
            
            ld=coolExp(0.55,1.06,mu_fibretype);
            ls=coolExp(77.5e-6*100,113e-6*100,mu_fibretype);
            
            rd=coolExp(41.5e-6*100,92.5e-6*100,mu_fibretype)/2;
            rs=coolExp(77.5e-6*100,113e-6*100,mu_fibretype)/2;
            
            c(1,:) = 2*pi*rd.*ld./Rmd;   % para.Gld
            c(2,:) = 4*2*pi*rs.*ls;      % para.Gkf
            c(3,:) = 16*2*pi*rs.*ls;     % para.Gks
            c(4,:) = 30*2*pi*rs.*ls;     % para.Gna
            c(5,:) = 2./(Ri*ld./(pi*rd.^2)+Ri*ls./(pi*rs.^2));     % para.Gc
            c(6,:) = 2*pi*rs.*ls./Rms;   % para.Gls
            c(7,:) = 2*pi*rd.*ld*Cm;     % para.Cd
            c(8,:) = 2*pi*rs.*ls*Cm;     % para.Cs
            s = ones(size(mu_fibretype));
            c(9,:) = 120*s;     % para.Vna
            c(10,:) = -10*s;     % para.Vk
            c(11,:) = 0*s;     % para.Vl
            
            function v = coolExp(a,b,mu)
                v = exp(log(100)*mu)*(b-a)/100 + a;
            end
        end
        
        dy = dydt_sarcomere(this, y, t);
        
        initSarcoConst(this);
        
        J = Jac_Sarco(this, y, t);
    end
end
