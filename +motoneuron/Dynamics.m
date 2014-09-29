classdef Dynamics < dscomponents.ACoreFun
    % FibreDynamics: Class for nonlinear dynamics of muscle fibre compound
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
   
    properties(Access=private)
        % mu-dep const
        mc;
    end
    
    methods
        function this = Dynamics(sys)
            this = this@dscomponents.ACoreFun(sys);
            this.initJSparsityPattern;
            this.fDim = 6;
            this.xDim = this.fDim;
        end
        
        function prepareSimulation(this, mu)
            if this.System.Model.FibreTypeDepMaxMeanCurrent
                mu(2) = min(polyval(this.System.upperlimit_poly,mu(1)),mu(2));
            end
            prepareSimulation@dscomponents.ACoreFun(this, mu);
            this.mc = this.getMotoConst(mu(1));
        end
        
        function dy = evaluate(this, y, ~, ~)
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
            
            c = this.mc;

            dy = zeros(size(y));

            % dendrites
            dy(1,:) = (-c(1,:).*(y(1,:)-c(11,:))-c(5,:).*(y(1,:)-y(2,:)))./c(7,:);
            
            % soma
            dy(2,:) = (-c(6,:).*(y(2,:)-c(11,:))-c(5,:).*(y(2,:)-y(1,:))...
                -c(4,:).*y(3,:).^3.*y(4,:).*(y(2,:)-c(9,:))...
                -c(2,:).*y(5,:).^4.*(y(2,:)-c(10,:))...
                -c(3,:).*y(6,:).^2.*(y(2,:)-c(10,:)))./c(8,:);
            
            % the four gating variables
            dy(3,:) = 0.32*(13-y(2,:))./(exp((13-y(2,:))/5)-1).*(1-y(3,:))...
                -0.28*(y(2,:)-40)./(exp((y(2,:)-40)/5)-1).*(y(3,:));
            dy(4,:) = 0.128*(exp((17-y(2,:))/18)).*(1-y(4,:))-4./(exp((40-y(2,:))/5)+1).*(y(4,:));
            dy(5,:) = 0.032*(15-y(2,:))./(exp((15-y(2,:))/5)-1).*(1-y(5,:))...
                -0.5*(exp((10-y(2,:))/40)).*(y(5,:));
            dy(6,:) = 3.5./(exp((55-y(2,:))/4)+1).*(1-y(6,:))-0.025*(y(6,:));
        end
        
        function dy = evaluateCoreFun(this, y, t, mu)%#ok
            error('evaluate is overridden directly.');
        end
        
        function J = getStateJacobian(this, y, ~, ~)
            % Evaluates the full Jacobian of the nonlinear core function.
            % Vectorized evaluation is not tested yet.
            %
            % Parameters:
            % y: @type matrix<double>, each column is one state vector
            % t: current time @type rowvec<double>
            % mu: fibre type parameter, 0 = slow twitch fibre, 1 = fast twitch fibre, @type rowvec<double>
            
            c = this.mc;
            J = sparse(6, 6);

            J(1,1:6:end) = -(c(1,:) + c(5,:))./c(7,:);
            J(2,1:6:end) = c(5,:)./c(8,:);
            J(1,2:6:end) = c(5,:)./c(7,:);
            J(2,2:6:end) = -(c(4,:).*y(4,:).*y(3,:)^3 + c(2,:).*y(5,:)^4 + c(3,:).*y(6,:)^2 + c(5,:) + c(6,:))./c(8,:);
            J(3,2:6:end) = (8.*(y(3,:) - 1))./(25.*(exp(13./5 - y(2,:)./5) - 1)) - (7.*y(3,:))./(25.*(exp(y(2,:)./5 - 8) - 1)) + (y(3,:).*exp(y(2,:)./5 - 8).*((7.*y(2,:))./25 - 56./5))./(5.*(exp(y(2,:)./5 - 8) - 1)^2) + (exp(13./5 - y(2,:)./5).*((8.*y(2,:))./25 - 104./25).*(y(3,:) - 1))./(5.*(exp(13./5 - y(2,:)./5) - 1)^2);
            J(4,2:6:end) = (8.*exp(17./18 - y(2,:)./18).*(y(4,:) - 1))./1125 - (4.*y(4,:).*exp(8 - y(2,:)./5))./(5.*(exp(8 - y(2,:)./5) + 1)^2);
            J(5,2:6:end) = (y(5,:).*exp(1./4 - y(2,:)./40))./80 + (4.*(y(5,:) - 1))./(125.*(exp(3 - y(2,:)./5) - 1)) + (exp(3 - y(2,:)./5).*((4.*y(2,:))./125 - 12./25).*(y(5,:) - 1))./(5.*(exp(3 - y(2,:)./5) - 1)^2);
            J(6,2:6:end) = -(7.*exp(55./4 - y(2,:)./4).*(y(6,:) - 1))./(8.*(exp(55./4 - y(2,:)./4) + 1)^2);
            J(2,3:6:end) = (3.*c(4,:).*y(3,:)^2.*y(4,:).*(c(9,:) - y(2,:)))./c(8,:);
            J(3,3:6:end) =  ((8.*y(2,:))./25 - 104./25)./(exp(13./5 - y(2,:)./5) - 1) - ((7.*y(2,:))./25 - 56./5)./(exp(y(2,:)./5 - 8) - 1);
            J(2,4:6:end) = (c(4,:).*y(3,:)^3.*(c(9,:) - y(2,:)))./c(8,:);
            J(4,4:6:end) =  - (16.*exp(17./18 - y(2,:)./18))./125 - 4./(exp(8 - y(2,:)./5) + 1);
            J(2,5:6:end) = (4.*c(2,:).*y(5,:)^3.*(c(10,:) - y(2,:)))./c(8,:);
            J(5,5:6:end) = ((4.*y(2,:))./125 - 12./25)./(exp(3 - y(2,:)./5) - 1) - exp(1./4 - y(2,:)./40)./2;
            J(2,6:6:end) = (2.*c(3,:).*y(6,:).*(c(10,:) - y(2,:)))./c(8,:);
            J(6,6:6:end) = - 7./(2.*(exp(55./4 - y(2,:)./4) + 1)) - 1./40;
        end
        
        function copy = clone(this)
            % Create new instance
            copy = models.motoneuron.Dynamics(this.System);
            
            % Call superclass clone (for deep copy)
            copy = clone@dscomponents.ACoreFun(this, copy);
            
            % Copy local properties
            copy.mc = this.mc;
        end
        
        function initJSparsityPattern(this)
            % Initializes the Sparsity pattern of the Jacobian
            i = [1,1,2,2,2,2,2,2,3,3,4,4,5,5,6,6];
            j = [1,2,1,2,3,4,5,6,2,3,2,4,2,5,2,6];
            this.JSparsityPattern = logical(sparse(i,j,true,6,6));
        end
    end
    
    methods(Access=private)
        function c = getMotoConst(~, moto_mu)
            % getMotoConst: private getter function for motoneuron constants
            %
            % moto_mu is assumed to be in [0,1]

            Cm=1;
            Ri=70/1000;
            c = zeros(11,length(moto_mu));
            Rmd = 14.4+6.05-coolExp(6.05,14.4,moto_mu);  % cf. Cisi and Kohn 2008, Table 2, page 7
            Rms=1.15+0.65-coolExp(0.65,1.15,moto_mu);

            ld=coolExp(0.55,1.06,moto_mu);
            ls=coolExp(77.5e-6*100,113e-6*100,moto_mu);

            rd=coolExp(41.5e-6*100,92.5e-6*100,moto_mu)/2;
            rs=coolExp(77.5e-6*100,113e-6*100,moto_mu)/2;

            c(1,:) = 2*pi*rd.*ld./Rmd;   % para.Gld
            c(2,:) = 4*2*pi*rs.*ls;      % para.Gkf
            c(3,:) = 16*2*pi*rs.*ls;     % para.Gks
            c(4,:) = 30*2*pi*rs.*ls;     % para.Gna
            c(5,:) = 2./(Ri*ld./(pi*rd.^2)+Ri*ls./(pi*rs.^2));     % para.Gc
            c(6,:) = 2*pi*rs.*ls./Rms;   % para.Gls
            c(7,:) = 2*pi*rd.*ld*Cm;     % para.Cd
            c(8,:) = 2*pi*rs.*ls*Cm;     % para.Cs
            s = ones(size(moto_mu));
            c(9,:) = 120*s;     % para.Vna
            c(10,:) = -10*s;     % para.Vk
            c(11,:) = 0*s;     % para.Vl

            function v = coolExp(a,b,mu)
                v = exp(log(100)*mu)*(b-a)/100 + a;
            end
        end
    end
end
