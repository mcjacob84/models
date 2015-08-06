classdef Force < dscomponents.ACoreFun
% Force: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2015-08-06
%
% @new{0,7,dw,2015-08-06} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        
    end
    
    methods
        function this = Force(sys)
            this = this@dscomponents.ACoreFun(sys);
            this.xDim = 4;
            this.fDim = 2;
            this.JSparsityPattern = sparse(true(2,4));
        end
        
        function dv = evaluate(this, x_v, t)           
            m = this.System.Model;
            mu = this.mu;
            
            gr = m.gradgreen(x_v(1:2));
            % Normale berechnen
            normale = [gr; -1]/norm([gr; -1]);
            % Normalenkraft
            Fn = normale * -mu(7)*normale(3);
            % Hangabtriebskraft
            Fb = [0;0;-mu(7)] - Fn;
            % Dreht den Fb-Vektror in die x,y-Richtung (kleiner Optimierungsgedanke)
            %Fb = -gr/norm(gr)*norm(Fb);
            
            %% MIT Reibung
            % Reibungskonstante * |F_n| * normierter Geschwindigkeitsvektor * zeiteinheit
            %Fr = mu*norm([gr;1]) * (v / norm(v)); %ALT
            Fr = mu(8)*norm(Fn) * (x_v(3:4) / norm(x_v(3:4)));
            % Gesamtkraft; wir lassen die z-Komponente der Hangabtriebskraft
            % ausser acht, da diese auf einem Putting-Green erwartungsgemäß
            % klein ist.
            dv = Fb([1 2])-Fr;
            
            %% OHNE Reibung
            %dv = Fb;
        end
        
        function evaluateCoreFun(~)
            error('evaluate is overridden directly.');
        end
    end
    
end