classdef MuscleFibreSystem < models.BaseDynSystem
% MuscleFibreSystem: The global dynamical system used within the MuscleFibreModel
%
% Contains the FibreDynamics, Linear Diffusion and InputConv components as well as a
% ConstInitialValue.
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
        % The number of sarcomer cells in this fibre model
        N;
    end
    
    properties (Constant)
        dm = 6; % Dimension of motoneuron part
        ds = 9; % Dimension of spindle part
        dsa = 58; % Dimension of single sarcomer cell part
    end
    
    methods
        function this = MuscleFibreSystem(model, N)
            % Call superclass constructor
            this = this@models.BaseDynSystem(model);
            % Set local variables
            this.N = N;
            
            %% Set system components
            % Core nonlinearity
            this.f = models.muscle.FibreDynamics(this);
            
            % Linear diffusion part
            this.A = dscomponents.LinearCoreFun(this.assembleA);
            
            % Linear input B for motoneuron
            this.B = dscomponents.LinearInputConv(this.assembleB);
            
            % Constant initial values
            this.x0 = dscomponents.ConstInitialValue(this.initStates);
        end
    end
    
    methods(Access=private)
        function A = assembleA(this)
            n0=this.dm+this.ds;  % size of upper left zeros block
            sar=this.dsa; % dimension of single sarcomer part
            size=this.dm+this.ds+sar*this.N;
            
            i=[n0+1+sar:sar:size,n0+1:sar:size,n0+1:sar:size-sar];    % row index
            j=[n0+1:sar:size-sar,n0+1:sar:size,n0+1+sar:sar:size];   % column index
            s=[ones(1,this.N-1),-1,-2*ones(1,this.N-2),-1,ones(1,this.N-1)];   % values
            A=sparse(i,j,s,size,size);
        end
        
        function B = assembleB(this)
            B=zeros(this.dm+this.ds+this.dsa*this.N,1);
            % B(2,:)=1/para.Cs
        end
        
        function x0 = initStates(this)
             % x0 = 0; noch n�tig?
             % return;
             x0 = [this.initNeuroStates;...
                this.initSpindleStates;...
                this.initSarcomerStates(this.N)];
        end
        
        function y = initNeuroStates(~)
            y = zeros(6,1);
        end
        
        function y = initSpindleStates(~)
            % copied from spindle_whole_2012_10_11.m, line 305 - 313
            y=zeros(9,1);
            y(1) = 0; %0.5*0.735;
            y(2) = 0; %0.5*0.735;
            y(3) = 0.000634066218078;
            y(4) = 0.000634066218078;
            y(5) = 0.000634066218078;
            y(6) = 0.020731324015575;
            y(7) = 0.020731324015575;
            y(8) = 0.020731324015575;
            y(9) = 0.95; %length
            %    y(:,10) = 0.0;  %velo
        end
        
        function y = initSarcomerStates(~,N)
            y=zeros(58,N);
            y(1,:) = -79.974;
            y(2,:) = -80.2;
            y(3,:) = 5.9;
            y(4,:) = 150.9;
            y(5,:) = 5.9;
            y(6,:) = 12.7;
            y(7,:) = 133;
            y(8,:) = 133;
            y(9,:) = 0.009466;
            y(10,:) = 0.9952;
            y(11,:) = 0.0358;
            y(12,:) = 0.4981;
            y(13,:) = 0.581;
            y(14,:) = 0.009466;
            y(15,:) = 0.9952;
            y(16,:) = 0.0358;
            y(17,:) = 0.4981;
            y(18,:) = 0.581;
            y(19,:) = 0;
            y(20,:) = 0;
            y(21,:) = 0;
            y(22,:) = 0;
            y(23,:) = 0;
            y(24,:) = 1;
            y(25,:) = 0;
            y(26,:) = 0;
            y(27,:) = 0;
            y(28,:) = 0;
            y(29,:) = 0.1;
            y(30,:) = 1500;
            y(31,:) = 0.1;
            y(32,:) = 1500;
            y(33,:) = 25;
            y(34,:) = 615;
            y(35,:) = 615;
            y(36,:) = 811;
            y(37,:) = 811;
            y(38,:) = 16900;
            y(39,:) = 16900;
            y(40,:) = 0.4;
            y(41,:) = 0.4;
            y(42,:) = 7200;
            y(43,:) = 7200;
            y(44,:) = 799.6;
            y(45,:) = 799.6;
            y(46,:) = 1000;
            y(47,:) = 1000;
            y(48,:) = 3;
            y(49,:) = 0.8;
            y(50,:) = 1.2;
            y(51,:) = 3;
            y(52,:) = 0.3;
            y(53,:) = 0.23;
            y(54,:) = 0.23;
            y(55,:) = 0.23;
            y(56,:) = 0.23;
            y(57,:) = 0.0;              % x_1 [micrometre]
            %    y(:,58) = CONSTANTS(:,106); % x_2 [micrometre]
            y(58,:) = 0.05;
            y=y(:);
        end
    end
    
end