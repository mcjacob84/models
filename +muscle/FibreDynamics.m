classdef FibreDynamics < dscomponents.ACompEvalCoreFun
    % FibreDynamics: Class for nonlinear dynamics of muscle fibre compound
    %
    % 
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
        system;
    end
    
    properties(SetAccess=private)
        MotoConst;
        spindleConst;
        SarcoConst;
    end
    
    methods
        function this = FibreDynamics(system)
            this.system = system;
            % Init neuro params -> in eigene initMotoConst private methode packen
            %             para = 0;
            %
            %             this.neuropara = para;
            
            this.initMotoConst;
            this.initSpindleConst;
            mu = 0;  % fast and slow twitch fibers
            this.initSarcoConst(mu);
        end
        
        function dy = evaluate(this, y, t, mu)
            dm = this.system.dm;
            ds = this.system.ds;
            dsa = this.system.dsa;
            N = this.system.N;
            
            % todo anpassen für matrix-wertiges y, also 58Nx1000 -> 58 x 1000N und zurück
            % SarcoRates=this.SarcomerRates(reshape(y(dm+ds+1:end,:),dsa,N)', t, mu);
            SarcoRates = this.SarcomerRates(reshape(y(dm+ds+1:end,:),dsa,[]), t, mu);
            SarcoRates = reshape(SarcoRates, dsa*N,[]);
            
            %% Inner dynamics (with partial nonlinear linking)
            dy = [this.NeuroRates(y(1:dm,:), t, mu);
                this.spindlerates(y(dm+1:dm+ds,:), t, mu, y(1:dm,:));
                SarcoRates];
            
            %% Link of motoneuron to sarcomer cell
            % Fix link moto-sarcomer for to middle cell
            n_link=round(N/2);
            dy(dm+ds+n_link*dsa+1,:) = dy(dm+ds+n_link*dsa+1,:) + 1.3*y(2,:)/this.SarcoConst(1);
        end
        
        function dy = evaluateCoreFun(this, y, t, mu)
            error('evaluate is overridden directly.');
            
        end
        
        function copy = clone(this)
            % Create new instance
            copy = models.muscle.MuscleCoreFun(this.system);
            
            % Call superclass clone (for deep copy)
            copy = clone@dscomponents.ACompEvalCoreFun(this, copy);
            
            % Copy local properties
            % No local properties are to be copied here, as so far everything is done in the
            % constructor.
        end
    end
    
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
             % `f_i(\vx) = f_i(\vx(ends(i-1){:}ends(i)));` @type rowvec<integer>
             % idx: The indices of `\vx`-entries in the global `\vx` vector w.r.t the `i`-th
             % point, e.g. `xglobal(i-1:i+1) = \vx(ends(i-1):ends(i))` @type rowvec<integer>
             % self: The positions in the `\vx` vector that correspond to the `i`-th output
             % dimension, if applicable (usually `f_i` depends on `x_i`, but not necessarily)
             % @type rowvec<integer>
             % x: A matrix `\vX` with the state space locations `\vx_i` in its columns @type
             % matrix<double>
             % t: The corresponding times `t_i` for each state `\vx_i` @type rowvec<double>
             % mu: The corresponding parameters `\mu_i` for each state `\vx_i`, as column matrix
             % @type matrix<double>
             %
             % Return values:
             % fx: A matrix with pts-many component function evaluations `f_i(\vx)` as rows and as
             % many columns as `\vX` had.
        end
        
        function dfx = evaluateComponentPartialDerivatives(this, pts, ends, idx, deriv, self, x, t, mu, dfxsel)
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
    end
    
    methods(Access=private)
        function dy = NeuroRates(this, y, t, mu)
            
            % adapted from combined_cell_model.m, line 738 - 745
            para = this.MotoConst;
            
            statesSize = size(y);
            dy=zeros(statesSize);
            if ( statesSize(1) == 1)
                y = y';
                dy=dy';
            end
            
%             if size(dy)(1) ~= this.dm   ?????
%                 error('y has wrong dimension')
%             end
            
            % dendrites
            dy(1,:) = (-para.Gld.*(y(1,:)-para.Vl)-para.Gc.*(y(1,:)-y(2,:)))./para.Cd;
            % soma
            dy(2,:) = (-para.Gls.*(y(2,:)-para.Vl)-para.Gc.*(y(2,:)-y(1,:))...
                       -para.Gna.*y(3,:).^3.*y(4,:).*(y(2,:)-para.Vna)...
                       -para.Gkf.*y(5,:).^4.*(y(2,:)-para.Vk)...
                       -para.Gks.*y(6,:).^2.*(y(2,:)-para.Vk))./para.Cs;
 %tomo: where is input???
            % the four gating variables
            dy(3,:) = 0.32*(13-y(2,:))./(exp((13-y(2,:))/5)-1).*(1-y(3,:))...
                      -0.28*(y(2,:)-40)./(exp((y(2,:)-40)/5)-1).*(y(3,:));
            dy(4,:) = 0.128*(exp((17-y(2,:))/18)).*(1-y(4,:))-4./(exp((40-y(2,:))/5)+1).*(y(4,:));
            dy(5,:) = 0.032*(15-y(2,:))./(exp((15-y(2,:))/5)-1).*(1-y(5,:))...
                      -0.5*(exp((10-y(2,:))/40)).*(y(5,:));
            dy(6,:) = 3.5/(exp((55-y(2,:))/4)+1).*(1-y(6,:))-0.025*(y(6,:));
        end
        
        function dy = spindlerates(this, y, t, mu, y_moto)
            % adapted from spindle_whole_2012_10_11.m, line 319 - 356
            
            % mit �bergeben? gamma sp�ter aus Zustand y bekannt, L aus
            % Mechanik?
            %tomo: TODO Berechnung aus Parametern 
            gamma_dyn=0; 
            gamma_sta=0;
            
            %tomo TODO length/velocity/acceleration mapping
            L_ddot=0; % acceleration
            L_dot=0;  % velocity
            
            if (L_dot>=0)
                C=1;
            else
                C=0.42;
            end
            
            statesSize = size(y);
            dy=zeros(statesSize);
            if ( statesSize(1) == 1)
                y = y';
                dy=dy';
            end
            
            c = this.spindleConst;
            
            %tomo
            %    dy(:,10) = 0; %acceleration
            %    dy(:,9) = y(:,10); %velocity
            dy(9,:) = L_dot; %velocity
            
            dy(1,:) = (gamma_dyn^c(:,22)./(gamma_dyn^c(:,22)+c(:,21).^c(:,22)) - y(1,:))./c(:,20); % f_dyn_b1
            dy(2,:) = (gamma_sta^c(:,52)./(gamma_sta^c(:,52)+c(:,51).^c(:,52)) - y(2,:))./c(:,20); % f_stat_b2
            f_sta_chain = gamma_sta^c(:,82)./(gamma_sta^c(:,82)+c(:,81).^c(:,82)); % f_stat_chain
            
            % bag1
            beta_bag1  = c(:,4)+c(:,5).*y(1,:); %Beta
            gamma_bag1 = c(:,7).*y(1,:); %Gamma
            dy(3,:) = c(:,1)./c(:,3).*(C.*beta_bag1.*sign(L_dot-y(3,:)./c(:,1))*abs(L_dot-y(3,:)./c(:,1)).^c(:,15).*(y(9,:) ...
              -c(:,17)-y(6,:)./c(:,1)-c(:,16))+c(:,2).*(y(9,:)-c(:,17)-y(6,:)./c(:,1)-c(:,18))+c(:,3).*L_ddot+gamma_bag1-y(6,:));
            dy(6,:) = y(3,:);
            
            % bag2
            beta_bag2  = c(:,34)+c(:,35).*y(1,:)+c(:,36).*y(2,:); %Beta
            gamma_bag2 = c(:,37).*y(1,:)+c(:,38).*y(2,:); %Gamma
            dy(4,:) = c(:,31)./c(:,33).*(C.*beta_bag2.*sign(L_dot-y(4,:)./c(:,31))*abs(L_dot-y(4,:)./c(:,31)).^c(:,45).*(y(9,:)- ...
              c(:,47)-y(7,:)./c(:,31)-c(:,46))+c(:,32).*(y(9,:)-c(:,47)-y(7,:)./c(:,31)-c(:,48))+c(:,33).*L_ddot+gamma_bag2-y(7,:));
            dy(7,:) = y(4,:);
            
            % chain
            beta_chain  = c(:,64)+c(:,65).*y(1,:)+c(:,66).*f_sta_chain ; %Beta
            gamma_chain = c(:,67).*y(1,:)+c(:,68).*f_sta_chain; %Gamma
            dy(5,:)  = c(:,61)./c(:,63).*(C.*beta_chain.*sign(L_dot-y(5,:)./c(:,61))*abs(L_dot-y(5,:)./c(:,61)).^c(:,75).*(y(9,:)- ...
              c(:,77)-y(8,:)./c(:,61)-c(:,76))+c(:,62).*(y(9,:)-c(:,77)-y(8,:)./c(:,61)-c(:,78))+c(:,63).*L_ddot+gamma_chain-y(8,:));
            dy(8,:)  = y(5,:);

        end
        
        function dy = SarcomerRates(this, y, t, mu)
            % beta=STATES(2) (link to Neuro Modell)
            
            statesSize = size(y);
            dy = zeros(statesSize);
            %alg = zeros(77,statesSize(2));
            
            statesColumnCount = statesSize(2);
            if ( statesColumnCount == 1)
                y = y';
                alg = zeros(1, 77);
                dy = dy';
            else
                statesRowCount = statesSize(1);
                alg = zeros(statesRowCount, 77);
            end
            
            %tomo
            a=0; % ursprünglich in initSarcoConsts() in Sarco Modell
                 % a = velocity, welche? ist die in Zustand y?
                 % used in line ~412
            
      %             if size(dy)(2) ~= this.dsa   ?????
      %                 error('y has wrong dimension')
      %             end
      
            c = this.SarcoConst;
            
            % ALGEBRAIC(:,33) = 1.3*STATES(:,60);
            % alg(:,33) = 1.3*beta; % neuro-sarco link now in evaluate

            % tomo TODO - when coupling to mechanics --> velocity mapping
            % L_S = L_S_0
            alg(:,73) = c(:,66);
            
            % T_tot = (a*L_S^4+b*L_S^3+c*L_S^2+d*L_S+e)*T_tot_0 --- a fourth order polynomial --- cf. mod_num_of_XB.m
            alg(:,74) = piecewise({alg(:,73)>=0.4&alg(:,73)<2.0,(0.073164682539683.*(2*alg(:,73)).^4-0.702380952380952.*(2*alg(:,73)).^3+1.950644841269842.*(2*alg(:,73)).^2-1.271666666666668.*(2*alg(:,73))+0.098571428571430).*c(:,70)},0.0);
            % g = g_0 * exp{ sigma * (x_2-x_0)^2 }
            alg(:,75) = c(:,93).*exp(c(:,110).*(y(:,58)-c(:,106)).^2);
            % f_new = f_o { 1 + lambda_A1 * [exp(x_1/x_0 * (nu-1)) - 1] + lambda_A2 * [exp(x_2/x_0*(nu-1)) - 1]}^2
            alg(:,77) = c(:,89) .* (1 + y(:,52)./alg(:,74) .* (exp(y(:,57)./c(:,106).*(c(:,109)-1))-1) + y(:,53)./alg(:,74) .* (exp(y(:,58)./c(:,106).*(c(:,109)-1))-1)).^2;
            
            dy(:,29) = (((( ( c(:,99).*(y(:,19)+y(:,20)+y(:,21)+y(:,22)+y(:,23))).*((y(:,30) - y(:,29))./c(:,102)) -  c(:,61).*((y(:,29)./(y(:,29)+c(:,62)))./c(:,102)))+ c(:,63).*((y(:,30) - y(:,29))./c(:,102)))+  - c(:,64).*((y(:,29) - y(:,31))./c(:,102)))+ - ( ( c(:,71).*y(:,29)).*((c(:,73)+ - y(:,34))+ - y(:,36))+  - c(:,72).*y(:,34)))+ - ( ( c(:,79).*y(:,29)).*y(:,44)+  - c(:,80).*y(:,40));
            dy(:,30) = (((  - ( c(:,99).*(y(:,19)+y(:,20)+y(:,21)+y(:,22)+y(:,23))).*((y(:,30) - y(:,29))./c(:,104))+ c(:,61).*((y(:,29)./(y(:,29)+c(:,62)))./c(:,104)))+  - c(:,63).*((y(:,30) - y(:,29))./c(:,104)))+  - c(:,65).*((y(:,30) - y(:,32))./c(:,104)))+ - ( ( c(:,76).*y(:,30)).*(c(:,78) - y(:,38))+  - c(:,77).*y(:,38));
            dy(:,32) = ((( c(:,61).*((y(:,31)./(y(:,31)+c(:,62)))./c(:,105))+  - c(:,63).*((y(:,32)+ ...
                - y(:,31))./c(:,105)))+ c(:,65).*((y(:,30)+ - y(:,32))./c(:,105)))+ - ( ( c(:,76).*y(:,32)).*(c(:,78)+ - y(:,39))+ ...
                - c(:,77).*y(:,39))) -  (1000.00./1.00000).*( c(:,96).*( y(:,55).*(0.00100000./1.00000).*y(:,32) - c(:,98)).*(piecewise({ y(:,55).*(0.00100000./1.00000).*y(:,32) - c(:,98)>0.00000, 1.00000 }, 0.00000)).*(0.00100000./1.00000).*y(:,55).*y(:,32) -  c(:,97).*y(:,56).*(c(:,98) -  y(:,55).*(0.00100000./1.00000).*y(:,32)).*(piecewise({c(:,98) -  y(:,55).*(0.00100000./1.00000).*y(:,32)>0.00000, 1.00000 }, 0.00000)));
            dy(:,34) =  ( c(:,71).*y(:,29)).*((c(:,73)+ - y(:,34))+ - y(:,36))+  - c(:,72).*y(:,34);
            dy(:,35) =  ( c(:,71).*y(:,31)).*((c(:,73)+ - y(:,35))+ - y(:,37))+  - c(:,72).*y(:,35);
            dy(:,36) =  ( c(:,74).*(c(:,73)+ - y(:,34)+ - y(:,36))).*y(:,46)+  - c(:,75).*y(:,36);
            dy(:,37) =  ( c(:,74).*(c(:,73)+ - y(:,35)+ - y(:,37))).*y(:,47)+  - c(:,75).*y(:,37);
            dy(:,38) =  ( c(:,76).*y(:,30)).*(c(:,78)+ - y(:,38))+  - c(:,77).*y(:,38);
            dy(:,39) =  ( c(:,76).*y(:,32)).*(c(:,78)+ - y(:,39))+  - c(:,77).*y(:,39);
            dy(:,40) = ( ( c(:,79).*y(:,29)).*y(:,44)+  - c(:,80).*y(:,40))+  - c(:,83).*((y(:,40)+ - y(:,41))./c(:,102));
            dy(:,41) = ( ( c(:,79).*y(:,31)).*y(:,45)+  - c(:,80).*y(:,41))+ c(:,83).*((y(:,40)+ - y(:,41))./c(:,103));
            dy(:,42) = ( ( c(:,81).*y(:,46)).*y(:,44)+  - c(:,82).*y(:,42))+  - c(:,83).*((y(:,42)+ - y(:,43))./c(:,102));
            dy(:,43) = ( ( c(:,81).*y(:,47)).*y(:,45)+  - c(:,82).*y(:,43))+ c(:,83).*((y(:,42)+ - y(:,43))./c(:,103));
            dy(:,44) = ( - ( ( c(:,79).*y(:,29)).*y(:,44)+  - c(:,80).*y(:,40))+ - ( ( c(:,81).*y(:,46)).*y(:,44)+  - c(:,82).*y(:,42)))+  - c(:,83).*((y(:,44)+ - y(:,45))./c(:,102));
            dy(:,45) = ( - ( ( c(:,79).*y(:,31)).*y(:,45)+  - c(:,80).*y(:,41))+ - ( ( c(:,81).*y(:,47)).*y(:,45)+  - c(:,82).*y(:,43)))+ c(:,83).*((y(:,44)+ - y(:,45))./c(:,103));
            dy(:,46) = ( - ( ( c(:,74).*(c(:,73)+ - y(:,34)+ - y(:,36))).*y(:,46)+  - c(:,75).*y(:,36))+ - ( ( c(:,81).*y(:,46)).*y(:,44)+  - c(:,82).*y(:,42)))+  - c(:,84).*((y(:,46)+ - y(:,47))./c(:,102));
            dy(:,47) = ( - ( ( c(:,74).*(c(:,73)+ - y(:,35)+ - y(:,37))).*y(:,47)+  - c(:,75).*y(:,37))+ - ( ( c(:,81).*y(:,47)).*y(:,45)+  - c(:,82).*y(:,43)))+ c(:,84).*((y(:,46)+ - y(:,47))./c(:,103));
            dy(:,48) = (( ( c(:,68).*y(:,31)).*y(:,33)+  - c(:,69).*y(:,48))+  - c(:,87).*y(:,48))+ c(:,88).*y(:,51);
            dy(:,50) = (((( c(:,68).*y(:,31).*y(:,49)+  - c(:,69).*y(:,50))+ c(:,85).*y(:,33))+  - c(:,86).*y(:,50))+ (  - c(:,68).*y(:,31)).*y(:,50))+ c(:,69).*y(:,51);
            dy(:,51) = ((((( c(:,68).*y(:,31).*y(:,50)+  - c(:,69).*y(:,51))+ c(:,87).*y(:,48))+  - c(:,88).*y(:,51))+  - alg(:,77).*y(:,51))+ c(:,90).*y(:,52))+ alg(:,75).*y(:,53);
            dy(:,52) = (( alg(:,77).*y(:,51)+  - c(:,90).*y(:,52))+ c(:,92).*y(:,53))+  - c(:,91).*y(:,52);
            dy(:,53) = (  - c(:,92).*y(:,53)+ c(:,91).*y(:,52))+  - alg(:,75).*y(:,53);
            dy(:,54) =  (0.00100000./1.00000).*( c(:,91).*y(:,52) -  c(:,92).*y(:,53))+ -1.00000.*c(:,94).*y(:,54)+ -1.00000.*c(:,95).*((y(:,54) - y(:,55))./c(:,103));
            dy(:,55) =  c(:,95).*((y(:,54) - y(:,55))./c(:,105)) -  1.00000.*( c(:,96).*( y(:,55).*(0.00100000./1.00000).*y(:,32) - c(:,98)).*(piecewise({ y(:,55).*(0.00100000./1.00000).*y(:,32) - c(:,98)>0.00000, 1.00000 }, 0.00000)).*(0.00100000./1.00000).*y(:,55).*y(:,32) -  c(:,97).*y(:,56).*(c(:,98) -  y(:,55).*(0.00100000./1.00000).*y(:,32)).*(piecewise({c(:,98) -  y(:,55).*(0.00100000./1.00000).*y(:,32)>0.00000, 1.00000 }, 0.00000)));
            dy(:,56) =  1.00000.*( c(:,96).*( y(:,55).*(0.00100000./1.00000).*y(:,32) - c(:,98)).*(piecewise({ y(:,55).*(0.00100000./1.00000).*y(:,32) - c(:,98)>0.00000, 1.00000 }, 0.00000)).*(0.00100000./1.00000).*y(:,55).*y(:,32) -  c(:,97).*y(:,56).*(c(:,98) -  y(:,55).*(0.00100000./1.00000).*y(:,32)).*(piecewise({c(:,98) -  y(:,55).*(0.00100000./1.00000).*y(:,32)>0.00000, 1.00000 }, 0.00000)));
            alg(:,13) = piecewise({alg(:,74)+ - y(:,33)+ - y(:,48)+ - y(:,49)+ - y(:,50)+ - y(:,51)+ - y(:,52)+ - y(:,53)>=0.00000, alg(:,74)+ - y(:,33)+ - y(:,48)+ - y(:,49)+ - y(:,50)+ - y(:,51)+ - y(:,52)+ - y(:,53) }, 0.00000);
            dy(:,31) = ((((  - c(:,61).*((y(:,31)./(y(:,31)+c(:,62)))./c(:,103))+ c(:,63).*((y(:,32)+ - y(:,31))./c(:,103)))+ c(:,64).*((y(:,29) - y(:,31))./c(:,103)))+ - ((((((( c(:,68).*y(:,31).*alg(:,13)+  - c(:,69).*y(:,33))+ c(:,68).*y(:,31).*y(:,33))+  - c(:,69).*y(:,48))+ c(:,68).*y(:,31).*y(:,49))+  - c(:,69).*y(:,50))+ c(:,68).*y(:,31).*y(:,50))+  - c(:,69).*y(:,51)))+ - ( ( c(:,71).*y(:,31)).*(c(:,73)+ - y(:,35)+ - y(:,37))+  - c(:,72).*y(:,35)))+ - ( ( c(:,79).*y(:,31)).*y(:,45)+  - c(:,80).*y(:,41));
            dy(:,33) = (((( ( c(:,68).*y(:,31)).*alg(:,13)+  - c(:,69).*y(:,33))+ (  - c(:,68).*y(:,31)).*y(:,33))+ c(:,69).*y(:,48))+  - c(:,85).*y(:,33))+ c(:,86).*y(:,50);
            dy(:,49) = (( (  - c(:,68).*y(:,31)).*y(:,49)+ c(:,69).*y(:,50))+ c(:,85).*alg(:,13))+  - c(:,86).*y(:,49);
            alg(:,2) =  c(:,17).*((y(:,1) - c(:,22))./(1.00000 - (exp(( - ((y(:,1) - c(:,22))./c(:,33)))))));
            alg(:,15) =  c(:,20).*(exp(( - ((y(:,1) - c(:,22))./c(:,35)))));
            dy(:,9) =  alg(:,2).*(1.00000 - y(:,9)) -  alg(:,15).*y(:,9);
            alg(:,3) = 1.00000./(1.00000+(exp(((y(:,1) - c(:,26))./c(:,29)))));
            alg(:,16) =  1000.00.*(exp(( - ((y(:,1)+40.0000)./25.7500))));
            dy(:,10) = (alg(:,3) - y(:,10))./alg(:,16);
            alg(:,5) =  c(:,16).*((y(:,1) - c(:,21))./(1.00000 - (exp(( - ((y(:,1) - c(:,21))./c(:,32)))))));
            alg(:,18) =  c(:,19).*(exp(( - ((y(:,1) - c(:,21))./c(:,34)))));
            dy(:,11) =  alg(:,5).*(1.00000 - y(:,11)) -  alg(:,18).*y(:,11);
            alg(:,4) =  c(:,15).*(exp(( - ((y(:,1) - c(:,23))./c(:,30)))));
            alg(:,17) = c(:,18)./(1.00000+(exp(( - ((y(:,1) - c(:,23))./c(:,31))))));
            dy(:,12) =  alg(:,4).*(1.00000 - y(:,12)) -  alg(:,17).*y(:,12);
            alg(:,6) = 1.00000./(1.00000+(exp(((y(:,1) - c(:,25))./c(:,28)))));
            alg(:,19) = 8571.00./(0.200000+ 5.65000.*(((y(:,1)+c(:,49))./100.000) .^ 2.00000));
            dy(:,13) = (alg(:,6) - y(:,13))./alg(:,19);
            alg(:,7) =  c(:,17).*((y(:,2) - c(:,22))./(1.00000 - (exp(( - ((y(:,2) - c(:,22))./c(:,33)))))));
            alg(:,20) =  c(:,20).*(exp(( - ((y(:,2) - c(:,22))./c(:,35)))));
            dy(:,14) =  alg(:,7).*(1.00000 - y(:,14)) -  alg(:,20).*y(:,14);
            alg(:,8) = 1.00000./(1.00000+(exp(((y(:,2) - c(:,26))./c(:,29)))));
            alg(:,21) =  1.00000.*(exp(( - ((y(:,2)+40.0000)./25.7500))));
            dy(:,15) = (alg(:,8) - y(:,15))./alg(:,21);
            alg(:,10) =  c(:,16).*((y(:,2) - c(:,21))./(1.00000 - (exp(( - ((y(:,2) - c(:,21))./c(:,32)))))));
            alg(:,23) =  c(:,19).*(exp(( - ((y(:,2) - c(:,21))./c(:,34)))));
            dy(:,16) =  alg(:,10).*(1.00000 - y(:,16)) -  alg(:,23).*y(:,16);
            alg(:,9) =  c(:,15).*(exp(( - ((y(:,2) - c(:,23))./c(:,30)))));
            alg(:,22) = c(:,18)./(1.00000+(exp(( - ((y(:,2) - c(:,23))./c(:,31))))));
            dy(:,17) =  alg(:,9).*(1.00000 - y(:,17)) -  alg(:,22).*y(:,17);
            alg(:,11) = 1.00000./(1.00000+(exp(((y(:,2) - c(:,25))./c(:,28)))));
            alg(:,24) = 8571.00./(0.200000+ 5.65000.*(((y(:,2)+c(:,49))./100.000) .^ 2.00000));
            dy(:,18) = (alg(:,11) - y(:,18))./alg(:,24);
            alg(:,12) =  0.500000.*c(:,58).*(exp(((y(:,2) - c(:,60))./( 8.00000.*c(:,59)))));
            alg(:,25) =  0.500000.*c(:,58).*(exp(((c(:,60) - y(:,2))./( 8.00000.*c(:,59)))));
            dy(:,24) =   - c(:,55).*y(:,24)+ c(:,56).*y(:,19)+ -4.00000.*alg(:,12).*y(:,24)+ alg(:,25).*y(:,25);
            dy(:,19) =  c(:,55).*y(:,24)+  - c(:,56).*y(:,19)+( -4.00000.*alg(:,12).*y(:,19))./c(:,57)+ c(:,57).*alg(:,25).*y(:,20);
            dy(:,25) =  4.00000.*alg(:,12).*y(:,24)+  - alg(:,25).*y(:,25)+(  - c(:,55).*y(:,25))./c(:,57)+ c(:,57).*c(:,56).*y(:,20)+ -3.00000.*alg(:,12).*y(:,25)+ 2.00000.*alg(:,25).*y(:,26);
            dy(:,20) = ( c(:,55).*y(:,25))./c(:,57)+  - c(:,56).*c(:,57).*y(:,20)+( 4.00000.*alg(:,12).*y(:,19))./c(:,57)+  - c(:,57).*alg(:,25).*y(:,20)+( -3.00000.*alg(:,12).*y(:,20))./c(:,57)+ 2.00000.*c(:,57).*alg(:,25).*y(:,21);
            dy(:,26) =  3.00000.*alg(:,12).*y(:,25)+ -2.00000.*alg(:,25).*y(:,26)+(  - c(:,55).*y(:,26))./(c(:,57) .^ 2.00000)+ (c(:,57) .^ 2.00000).*c(:,56).*y(:,21)+ -2.00000.*alg(:,12).*y(:,26)+ 3.00000.*alg(:,25).*y(:,27);
            dy(:,21) = ( 3.00000.*alg(:,12).*y(:,20))./c(:,57)+ -2.00000.*c(:,57).*alg(:,25).*y(:,21)+( c(:,55).*y(:,26))./(c(:,57) .^ 2.00000)+  - c(:,56).*(c(:,57) .^ 2.00000).*y(:,21)+( -2.00000.*alg(:,12).*y(:,21))./c(:,57)+ 3.00000.*c(:,57).*alg(:,25).*y(:,22);
            dy(:,27) =  2.00000.*alg(:,12).*y(:,26)+ -3.00000.*alg(:,25).*y(:,27)+(  - c(:,55).*y(:,27))./(c(:,57) .^ 3.00000)+ c(:,56).*(c(:,57) .^ 3.00000).*y(:,22)+  - alg(:,12).*y(:,27)+ 4.00000.*alg(:,25).*y(:,28);
            dy(:,22) = ( c(:,55).*y(:,27))./(c(:,57) .^ 3.00000)+  - c(:,56).*(c(:,57) .^ 3.00000).*y(:,22)+( 2.00000.*alg(:,12).*y(:,21))./c(:,57)+ -3.00000.*alg(:,25).*c(:,57).*y(:,22)+(  - alg(:,12).*y(:,22))./c(:,57)+ 4.00000.*c(:,57).*alg(:,25).*y(:,23);
            dy(:,28) =  alg(:,12).*y(:,27)+ -4.00000.*alg(:,25).*y(:,28)+(  - c(:,55).*y(:,28))./(c(:,57) .^ 4.00000)+ c(:,56).*(c(:,57) .^ 4.00000).*y(:,23);
            dy(:,23) = ( alg(:,12).*y(:,22))./c(:,57)+ -4.00000.*c(:,57).*alg(:,25).*y(:,23)+( c(:,55).*y(:,28))./(c(:,57) .^ 4.00000)+  - c(:,56).*(c(:,57) .^ 4.00000).*y(:,23);
            alg(:,31) =  y(:,1).*((y(:,4) -  y(:,5).*(exp((( -1.00000.*c(:,7).*y(:,1))./( c(:,36).*c(:,37))))))./(1.00000 - (exp((( -1.00000.*c(:,7).*y(:,1))./( c(:,36).*c(:,37)))))));
            alg(:,14) =  (( c(:,36).*c(:,37))./c(:,7)).*(log((y(:,5)./y(:,4))));
            alg(:,38) =  y(:,5).*(exp(( (  - c(:,42).*alg(:,14)).*(c(:,7)./( c(:,36).*c(:,37))))));
            alg(:,39) =  c(:,41).*((alg(:,38) .^ 2.00000)./(c(:,43)+(alg(:,38) .^ 2.00000)));
            alg(:,40) = 1.00000 - ((1.00000+( c(:,44).*(1.00000+(alg(:,38) .^ 2.00000)./c(:,43)))./( (c(:,47) .^ 2.00000).*(exp((( 2.00000.*(1.00000 - c(:,42)).*y(:,1).*c(:,7))./( c(:,36).*c(:,37))))))) .^ -1.00000);
            alg(:,41) =  alg(:,39).*alg(:,40);
            alg(:,42) =  alg(:,41).*(piecewise({alg(:,31)>0.00000, 1.00000 }, 0.00000)).*(alg(:,31)./50.0000);
            alg(:,43) =  ( c(:,39).*(y(:,9) .^ 4.00000)).*y(:,10);
            alg(:,44) =  alg(:,43).*(alg(:,31)./50.0000);
            alg(:,48) =  (1.00000./7.00000).*((exp((y(:,8)./67.3000))) - 1.00000);
            alg(:,49) = (1.00000+ 0.120000.*(exp(( -0.100000.*y(:,1).*(c(:,7)./( c(:,36).*c(:,37))))))+ 0.0400000.*alg(:,48).*(exp(( - ( y(:,1).*(c(:,7)./( c(:,36).*c(:,37)))))))) .^ -1.00000;
            alg(:,50) =  c(:,7).*(c(:,48)./( ((1.00000+c(:,45)./y(:,5)) .^ 2.00000).*((1.00000+c(:,46)./y(:,6)) .^ 3.00000)));
            alg(:,51) =  alg(:,50).*alg(:,49);
            dy(:,5) = (alg(:,42)+alg(:,44)+c(:,13)+  - 2.00000.*alg(:,51))./( (1000.00./1.00000).*c(:,7).*c(:,6))+(y(:,3) - y(:,5))./c(:,11);
            alg(:,46) =  ( ( c(:,40).*(y(:,11) .^ 3.00000)).*y(:,12)).*y(:,13);
            alg(:,45) =  y(:,1).*((y(:,6) -  y(:,8).*(exp((( -1.00000.*c(:,7).*y(:,1))./( c(:,36).*c(:,37))))))./(1.00000 - (exp((( -1.00000.*c(:,7).*y(:,1))./( c(:,36).*c(:,37)))))));
            alg(:,47) =  alg(:,46).*(alg(:,45)./75.0000);
            dy(:,8) = (alg(:,47)+c(:,14)+ 3.00000.*alg(:,51))./( (1000.00./1.00000).*c(:,7).*c(:,6))+(y(:,7) - y(:,8))./c(:,12);
            alg(:,1) =  (1000.00./1.00000).*((y(:,1) - y(:,2))./c(:,3));
            alg(:,27) = 156.500./(5.00000+(exp(((  - c(:,7).*alg(:,14))./( c(:,36).*c(:,37))))));
            alg(:,28) = 156.500 -  5.00000.*alg(:,27);
            alg(:,35) =  y(:,1).*((alg(:,27) -  alg(:,28).*(exp((( c(:,7).*y(:,1))./( c(:,36).*c(:,37))))))./(1.00000 - (exp((( c(:,7).*y(:,1))./( c(:,36).*c(:,37)))))));
            alg(:,34) = 1.00000./(1.00000+(exp(((y(:,1) - c(:,24))./c(:,27)))));
            alg(:,36) =  c(:,38).*(alg(:,34) .^ 4.00000);
            alg(:,37) =  alg(:,36).*(alg(:,35)./45.0000);
            
            alg(:,52) = alg(:,37)+alg(:,42)+alg(:,44)+alg(:,47)+alg(:,51);  % + - alg(:,33);  % neuro-sarco link now in evaluate 
            dy(:,1) =  - ((alg(:,52)+alg(:,1))./c(:,1));
            alg(:,32) =  y(:,2).*((y(:,4) -  y(:,3).*(exp((( -1.00000.*c(:,7).*y(:,2))./( c(:,36).*c(:,37))))))./(1.00000 - (exp((( -1.00000.*c(:,7).*y(:,2))./( c(:,36).*c(:,37)))))));
            alg(:,26) =  (( c(:,36).*c(:,37))./c(:,7)).*(log((y(:,3)./y(:,4))));
            alg(:,57) =  y(:,3).*(exp(( (  - c(:,42).*alg(:,26)).*(c(:,7)./( c(:,36).*c(:,37))))));
            alg(:,58) =  c(:,41).*((alg(:,57) .^ 2.00000)./(c(:,43)+(alg(:,57) .^ 2.00000)));
            alg(:,59) = 1.00000 - ((1.00000+( c(:,44).*(1.00000+(alg(:,57) .^ 2.00000)./c(:,43)))./( (c(:,47) .^ 2.00000).*(exp((( 2.00000.*(1.00000 - c(:,42)).*y(:,2).*c(:,7))./( c(:,36).*c(:,37))))))) .^ -1.00000);
            alg(:,60) =  alg(:,58).*alg(:,59);
            alg(:,61) =  c(:,51).*alg(:,60).*(alg(:,32)./50.0000);
            alg(:,62) =  ( c(:,39).*(y(:,14) .^ 4.00000)).*y(:,15);
            alg(:,63) =  c(:,52).*alg(:,62).*(alg(:,32)./50.0000);
            alg(:,67) =  (1.00000./7.00000).*((exp((y(:,7)./67.3000))) - 1.00000);
            alg(:,68) = (1.00000+ 0.120000.*(exp(( -0.100000.*y(:,2).*(c(:,7)./( c(:,36).*c(:,37))))))+ 0.0400000.*alg(:,67).*(exp(( - ( y(:,2).*(c(:,7)./( c(:,36).*c(:,37)))))))) .^ -1.00000;
            alg(:,69) =  c(:,7).*(c(:,48)./( ((1.00000+c(:,45)./y(:,3)) .^ 2.00000).*((1.00000+c(:,46)./y(:,6)) .^ 3.00000)));
            alg(:,70) =  c(:,54).*alg(:,69).*alg(:,68);
            dy(:,4) =   - c(:,10).*((alg(:,61)+alg(:,63)+c(:,13)+  - 2.00000.*alg(:,70))./( (1000.00./1.00000).*c(:,7).*c(:,4))) - (alg(:,42)+alg(:,44)+c(:,13)+ -2.00000.*alg(:,51))./( (1000.00./1.00000).*c(:,7).*c(:,5));
            dy(:,3) = (alg(:,61)+alg(:,63)+c(:,13)+  - 2.00000.*alg(:,70))./( (1000.00./1.00000).*c(:,7).*c(:,4)) - (y(:,3) - y(:,5))./c(:,8);
            alg(:,65) =  ( ( c(:,40).*(y(:,16) .^ 3.00000)).*y(:,17)).*y(:,18);
            alg(:,64) =  y(:,2).*((y(:,6) -  y(:,7).*(exp((( -1.00000.*c(:,7).*y(:,2))./( c(:,36).*c(:,37))))))./(1.00000 - (exp((( -1.00000.*c(:,7).*y(:,2))./( c(:,36).*c(:,37)))))));
            alg(:,66) =  c(:,53).*alg(:,65).*(alg(:,64)./75.0000);
            dy(:,6) =   - c(:,10).*((alg(:,66)+c(:,14)+ 3.00000.*alg(:,70))./( (1000.00./1.00000).*c(:,7).*c(:,4))) - (alg(:,47)+c(:,14)+ 3.00000.*alg(:,51))./( (1000.00./1.00000).*c(:,7).*c(:,5));
            dy(:,7) = (alg(:,66)+c(:,14)+ 3.00000.*alg(:,70))./( (1000.00./1.00000).*c(:,7).*c(:,4)) - (y(:,7) - y(:,8))./c(:,9);
            alg(:,29) = 156.500./(5.00000+(exp(((  - c(:,7).*alg(:,26))./( c(:,36).*c(:,37))))));
            alg(:,30) = 156.500 -  5.00000.*alg(:,29);
            alg(:,54) =  y(:,2).*((alg(:,29) -  alg(:,30).*(exp((( c(:,7).*y(:,2))./( c(:,36).*c(:,37))))))./(1.00000 - (exp((( c(:,7).*y(:,2))./( c(:,36).*c(:,37)))))));
            alg(:,53) = 1.00000./(1.00000+(exp(((y(:,2) - c(:,24))./c(:,27)))));
            alg(:,55) =  c(:,38).*(alg(:,53) .^ 4.00000);
            alg(:,56) =  c(:,50).*alg(:,55).*(alg(:,54)./45.0000);
            alg(:,71) = alg(:,56)+alg(:,61)+alg(:,63)+alg(:,66)+alg(:,70);
            dy(:,2) =  - ((alg(:,71) - alg(:,1)./c(:,2))./c(:,1));
            % d/dt(x_1) = -(f_0*D_2/A_1 + h_p*A_2/A_1)*x_1 + h_p*A_2/A_1*(x_2-x_0) + velo/2
            dy(:,57) = -(alg(:,77).*y(:,51)./y(:,52)+c(:,92).*y(:,53)./y(:,52)).*y(:,57) + c(:,92).*y(:,53)./y(:,52).*(y(:,58)-c(:,106));% replaced by a + c(:,108)./2;
            % d/dt(x_2) = -h_0*A_1/A_2*(x_2-(x_1+x_0))  + velo/2
            dy(:,58) = -c(:,91).*y(:,52)./y(:,53).*(y(:,58)-(y(:,57)+c(:,106)));% replaces by a  + c(:,108)./2;

            %tomo USELESS HERE - ALG(:,72) is the active stress/force (we are interested in this component) - compute outside
            % Force = eta*(A_1/T_tot_0*x_1 + A_2/T_tot_0*x_2)
            %alg(:,72) = c(:,107).*(y(:,52)./c(:,70).*y(:,57)+y(:,53)./c(:,70).*y(:,58));

            %tomo TODO - for coupling to mechanics, the velocity has to be included
            % I think you can not have two equations with the same dy(i) on the left hand side !?!
            %dy(:,57) = dy(:,57) - a/2;
            %dy(:,58) = dy(:,58) - a/2;
            
            dy = dy';
            
            function x = piecewise(cases, default)
                set = 0;
                for i = 1:2:length(cases)
                    if (length(cases{i+1}) == 1)
                        x(cases{i} & ~set,:) = cases{i+1};
                    else
                        x(cases{i} & ~set,:) = cases{i+1}(cases{i} & ~set);
                    end
                    set = set | cases{i};
                    if(set), break, end
                end
                if (length(default) == 1)
                    x(~set,:) = default;
                else
                    x(~set,:) = default(~set);
                end
            end
        end
    end
    
    methods%(Access=private)
        
        function para = initMotoConst(this) %, mu)
            N = this.system.N_neuro;
            mu = 1:N;
            Cm=1;
            Ri=70/1000;
            Rmd = 14.4+6.05-coolExp(6.05,14.4,mu);  % cf. Cisi and Kohn 2008, Table 2, page 7
            Rms=1.15+0.65-coolExp(0.65,1.15,mu);
            
            ld=coolExp(0.55,1.06,mu);
            ls=coolExp(77.5e-6*100,113e-6*100,mu);

            rd=coolExp(41.5e-6*100,92.5e-6*100,mu)/2; 
            rs=coolExp(77.5e-6*100,113e-6*100,mu)/2;
            
            para.Gld=2*pi*rd.*ld./Rmd;
            para.Gkf=4*2*pi*rs.*ls;
            para.Gks=16*2*pi*rs.*ls;
            para.Gna=30*2*pi*rs.*ls;
            para.Gc=2./(Ri*ld./(pi*rd.^2)+Ri*ls./(pi*rs.^2));
            para.Gls=2*pi*rs.*ls./Rms;
            para.Cd=2*pi*rd.*ld*Cm;
            para.Cs=2*pi*rs.*ls*Cm;
            s = ones(size(mu));
            para.Vna=120*s;
            para.Vk=-10*s;
            para.Vl=0*s;
            
            this.MotoConst = para;
            
%             function vector=expdistribution(LL,UP,N)
%                 RR=UP-LL;
%                 a = log(100)/N;
%                 vector = exp(a*(1:N))/100*RR+LL;
%             end
            
%             function v = coolExp(a,b,mu)
%                 v = exp(log(100)*mu)*b/100 + a;     % fehlt da ein
%                 /NTOTAL in exp()?
%             end
            function v = coolExp(a,b,mu)
                v = exp(log(100)*mu/N)*b/100 + a;
            end
        
        end
        
        function initSpindleConst(this)
            
            % copied from spindle_whole_2012_10_11.m, line 217 - 291
            c = zeros(1,83); %y = [];
            
            %Non-existant for Bag_1: c(6), c(8), c(11), c(13), c(19)
            % c(9) replaced with C as an input parameter
            c(1:11) = [10.4649, 0.15, 0.0002, 0.0605, 0.2592, 0, 0.0289, 0, 0, 0, 0];
            c(12:22) = [0.0423, 0, 20000, 0.3, 0.46, 0.04, 0.76, 0, 0.149, 60, 2];
            
            %Non-existant for Bag_2: c(35), c(37); c(39) replaced with C as an input parameter 
            c(31:42) = [10.4649, 0.15, 0.0002, 0.0822, 0, -0.046, 0, 0.0636, 0, 0, 0.7, 0.0423];
            c(43:53) = [0.89, 10000, 0.3, 0.46, 0.04, 0.76, 0.04, 0.205, 60, 2, 7250];
            
            %Non-existant for Chain: c(65), c(67),c(80); c(69) replaced with C as an input parameter            
            c(61:72) = [10.4649, 0.15, 0.0002, 0.0822, 0, -0.069, 0, 0.0954, 0, 0, 0.7, 0.0423];
            c(73:83) = [0.89, 10000, 0.3, 0.46, 0.04, 0.76, 0.04, 0, 90, 2, 7250];
            
            this.spindleConst=c;
        end

        function initSarcoConst(this, mu)
            %slow twitch fibres
            c1 = zeros(1,110);
            c1(1:10)  = [0.58, 2.79, 150, 0.000001, 0.0025, 0.0005, 96485, 559, 559, 0.00174];
            c1(11:20) = [40229.885, 40229.885, 0.34, -0.43, 0.0081, 0.288, 0.0131, 4.38, 1.38, 0.067];
            c1(21:30) = [-46, -40, -45, 70, -68, -40, 150, 7.1, 7.5, 14.7];
            c1(31:40) = [9, 10, 7, 18, 40, 8314.41, 293, 3.275, 10.8, 134];
            c1(41:50) = [1.85, 0.4, 950, 1, 1, 13, 10, 0.0001656, 70, 0.1];
            c1(51:60) = [1.0, 0.45, 0.1, 0.1, 0.002, 1000, 0.2, 0.2, 4.5, -20];
            c1(61:70) = [2.4375, 1, 0.00004, 0.75, 0.75, 1.1, 0.5, 0.0885, 0.115, 50]; %tomo 2012/11/29: changed C(:,70) from 140 to 50
            c1(71:80) = [0, 0, 1500, 0, 0, 0.000004, 0.005, 31000, 0.15, 30];
            c1(81:90) = [0.0015, 0.15, 0.375, 1.5, 0, 0.15, 0.15, 0.05, 0.5, 5];    %    CONSTANTS(:,89) = 3.0; % tomo ????
            c1(91:100) = [0.08, 0.06, 0.04, 0.00000394, 0.00000362, 1, 0.0001, 6, 60, 0.950000.*c1(66).* pi.*(c1(67).^2.00000)];
            c1(101:105) = [0.0500000.*c1(66).* pi.*(c1(67).^2.00000), 0.0100000.*c1(100), 0.990000.*c1(100), 0.0100000.*c1(101), 0.990000.*c1(101)];
            c1(106) = 0.05;     % x_0 [micrometre]
            c1(107) = 1.E4;     % eta [micronewton_per_micrometre] tomo: NOT USED
            %tomo TODO - velocity mapping --> coupling to mechanics
            c1(108) = 0; %  old: c(108) = -a;       %-0.001;      % velo [micrometre_per_millisecond] %-0.001
            c1(109) = 3.4;  % nu >= 1.   nu = 1 --> no cooperative effects
            c1(110) = 1000; % sigma

            %fast twitch fibres
            c2 = zeros(1,110);
            c2(1:10)  = [1., 4.8, 150, 0.000001, 0.0025, 0.0005, 96485, 350, 350, 0.0032];
            c2(11:20) = [21875., 21875., 1.02, -1.29, 0.0081, 0.288, 0.0131, 4.38, 1.38, 0.067];
            c2(21:30) = [-46, -40, -45, 70, -78, -40, 150, 5.8, 7.5, 14.7];
            c2(31:40) = [9, 10, 7, 18, 40, 8314.41, 293, 19.65, 64.8, 804];
            c2(41:50) = [11.1, 0.4, 950, 1, 1, 13, 10, 0.000621, 90, 0.1];
            c2(51:60) = [1.0, 0.45, 0.1, 0.1, 0.002, 1000, 0.2, 0.2, 4.5, -20];
            c2(61:70) = [4.875, 1, 0.00002, 0.75, 0.75, 1.1, 0.5, 0.04425, 0.115, 140];
            c2(71:80) = [0.0417, 0.0005, 1500, 0.000033, 0.003, 0.000004, 0.005, 31000, 0.15, 30];
            
            c2(81:90) = [0.0015, 0.15, 0.375, 1.5, 0, 0.15, 0.15, 0.05, 1.5, 15];    %    CONSTANTS(:,89) = 3.0; % tomo ????
            c2(91:100) = [0.24, 0.18, 0.12, 0.00002867, 0.00000362, 1, 0.0001, 6, 300, 0.950000.*c2(66).* pi.*(c2(67) .^ 2.00000)];
            c2(101:105) = [0.0500000.*c2(66).* pi.*(c2(67) .^ 2.00000), 0.0100000.*c2(100), 0.990000.*c2(100), 0.0100000.*c2(101), 0.990000.*c2(101)];
            c2(106) = 0.05;     % x_0 [micrometre]
            c2(107) = 1.E4;     % eta [micronewton_per_micrometre] tomo: NOT USED
            %tomo
            c2(108) = 0; %  old: c(108) = -a;       %-0.001;      % velo [micrometre_per_millisecond] %-0.001
            c2(109) = 3.4;  % nu >= 1.   nu = 1 --> no cooperative effects
            c2(110) = 1000; % sigma

            %linear interpolation - small MUs are slow twtich type, large MUs are fast twitch type
            c = (1-mu)*c1 + mu*c2;
            
            this.SarcoConst = c;
        end
    end
    
end
