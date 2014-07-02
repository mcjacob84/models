classdef System < models.BaseDynSystem
% MotoSystem: The global dynamical system used within the MotoModel
%
%
% @author Daniel Wirtz @date 2012-11-22
%
% @new{0,7,dw,2013-07-05} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    properties(Transient)
        % debug output
        d = [];
    end
    
    properties (Constant)        
        dm = 6; % Dimension of motoneuron part
        
        % Membrane capacitance. Different values for different fibre types, 
        % due to different action potential propagation speeds
        % C_m is computed parameter-dependent.
        % These constants are for both slow and fast muscle models and are also used in the
        % first entry of the sarcomere constants computed in
        % models.muscle.FibreDynamics.initSarcoConst @type double
        C_m_slow = 0.58;
        C_m_fast = 1;
    end
    
    properties(SetAccess=private)
        % The NoiseGenerator to obtain the inputs u(t)      
        noiseGen;
    end
    
    methods
        function this = System(model)
            % Call superclass constructor
            this = this@models.BaseDynSystem(model);
            
            this.addParam('fibre_type', [0 1], 10);
            this.addParam('mean_current_factor', [.5 10], 10);
            
            % Setup noise input
            ng = models.motoneuron.NoiseGenerator;
            this.noiseGen = ng;
            this.Inputs{1} = @ng.getInput;
            
            %% Set system components
            % Core nonlinearity
            this.f = models.motoneuron.Dynamics(this);
            
            % Linear input B for motoneuron
            % input conversion matrix, depends on fibre type. Has only one
            % entry in second row.
            B = dscomponents.AffLinInputConv;
            % Base noise input mapping
            B.addMatrix('1./(pi*(exp(log(100)*mu(1,:))*3.55e-05 + 77.5e-4).^2)',...
                sparse(2,1,1,6,2));
            % Independent noise input mapping with µ_2 as mean current factor
            B.addMatrix('mu(2,:)./(pi*(exp(log(100)*mu(1,:))*3.55e-05 + 77.5e-4).^2)',...
                sparse(2,2,1,6,2));
            this.B = B;
            
            % Constant initial values
            this.x0 = dscomponents.ConstInitialValue(zeros(6,1));
        end
        
        function setConfig(this, mu, inputidx)
            setConfig@models.BaseDynSystem(this, mu, inputidx);
            this.noiseGen.setMu(mu);
            this.MaxTimestep = this.Model.dt*1000;            
        end
        
        function pm = plot(this, t, y, pm)
            if nargin < 4
                n = 2;
                if ~isempty(this.mu)
                    n = 3;
                end
                pm = PlotManager(false,1,n);
                pm.LeaveOpen = true;
            end
            
            % Input
            h = pm.nextPlot('input','Input currents (base + indep)','time','value');
            plot(h,t,this.noiseGen.getInput(t));
            
            % Effective input
            if ~isempty(this.mu)
                h = pm.nextPlot('eff_input','Effective input current','time','value');
                B = this.B.compose([],this.mu);
                plot(h,t,B(2,:)*this.noiseGen.getInput(t));
            end
            
            % Firing rate
            h = pm.nextPlot('moto','Motoneuron V_s','time','value');
            plot(h,t,y(2,:));
            
            if nargin < 4
                pm.done;
            end
        end
    end    
end
