classdef DLTNonlinearCoreFun < models.beam.DLTBaseCoreFun
% DLTNonlinearCoreFun: 
%
%
%
% @author Daniel Wirtz @date 2011-12-07
%
% @new{0,6,dw,2011-12-07} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(Access=private)
        R_big;
        currentx = [];
        B_big;
        B_big_const;
    end
    
    methods
        function this = DLTNonlinearCoreFun(system)
            this = this@models.beam.DLTBaseCoreFun(system);
            this.initialize;
        end
        
        function initialize(this)
            initialize@models.beam.DLTBaseCoreFun(this);
            % Initialize the constant part of B
            s = length(this.sys.Model.free);
            null = sparse(s,s);
            this.B_big_const = [null -speye(s); null null];
        end
        
        function fx = evaluateCoreFun(this, x, t, mu)%#ok
            % Daten nicht aktuell? (check inside)
            this.updateB(x, mu);
            % Funktion f auswerten
            fx = - this.B_big*x - this.R_big;
        end
        
        function J = getStateJacobian(this, x, t, mu)%#ok
            % Daten nicht aktuell? (check inside)
            this.updateB(x, mu);
            % Jacobi-Matrix ausgeben
            J = -this.B_big;
        end
        
        function pr = project(this, V, W)%#ok
            
        end
    end
    
    methods(Access=private)
        function updateB(this, x, mu)
            m = this.sys.Model;
            % Updates the x, J and R big versions (if needed)
            if isempty(this.currentx) || (norm(this.currentx - x) > 1e-8)
                [R_F, K_F] = this.weak_form(x);

                R = R_F(m.free);
                K = K_F(m.free, m.free);
                
                s = length(m.free);
                null = sparse(s,s);
                % Use K here to have even nonlinear C!
                K0 = this.sys.K0(m.free,m.free);
                C = mu(1)*K0 + mu(2)*this.sys.M_small;
                this.B_big = this.B_big_const + [null, null;
                                                 K, C];
                this.R_big = [0*R; R];
                this.currentx = x;
            end
        end
        
        function [R, K] = weak_form(this, u)
            % Nichtlineare Funktion (Schwache Form, R-f_s) und ihre Ableitung (nach den Freiheitsgraden) (tangentielle Steifigkeitsmatrix, K)
            m = this.sys.Model;
            data = m.data;
            
            % Steifigkeitsmatrix für u und T
            K = sparse(7 * data.num_knots, 7 * data.num_knots);
            % Residuumsvektor
            R = sparse(7 * data.num_knots, 1);

            % Tangentiale Steifigkeitsmatrix aufstellen und schwache Form mit der aktuellen Verschiebung auswerten
            el = m.Elements;
%             i = []; j = []; K = []; R = []; ir = [];
            for k=1:length(el)
                index_glob = el{k}.getGlobalIndices;
                [K_lok, R_lok] = el{k}.getLocalTangentials(u(index_glob));
                
                % Unklar welche alternative schneller ist; muss man in
                % späteren tests mal profilen. Idealerweise müssen die i,j
                % etc vektoren auch vorinitialisiert werden..
%                 [li,lj] = meshgrid(index_glob);
%                 i = [i; li(:)];
%                 j = [j; lj(:)];
%                 K = [K; K_lok(:)];
%                 
%                 ir = [ir; index_glob'];
%                 R = [R; R_lok];
                K(index_glob, index_glob) = K(index_glob, index_glob) + K_lok;
                R(index_glob) = R(index_glob) + R_lok;
            end
%             % Steifigkeitsmatrix für u und T
%             K = sparse(i,j,K,7 * data.num_knots, 7 * data.num_knots);
%             % Residuumsvektor
%             R = sparse(ir,ones(size(ir)),R, 7 * data.num_knots, 1);
        end
    end
end