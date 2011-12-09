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
        x_big = [];
        J_big;
    end
    
    methods
        function this = DLTNonlinearCoreFun(system)
            this = this@models.beam.DLTBaseCoreFun(system);
        end
        
        function fx = evaluateCoreFun(this, x, t, mu)%#ok
            % Daten nicht aktuell? (check inside)
            this.updatexJ(x);
            % Funktion f auswerten
            fx = this.f_big - this.B_big*x - this.R_big;
        end
        
        function J = getStateJacobian(this, x, t, mu)%#ok
            % Daten nicht aktuell? (check inside)
            this.updatexJ(x);
            % Jacobi-Matrix ausgeben
            J = this.J_big;
        end
        
        function pr = project(this, V, W)%#ok
            
        end
        
%         function prepareConstants(this, ~)
%             prepareConstants@models.beam.DLTBaseCoreFun(this);
%             
%             % Assemble constant part of B_big matrix (varying K for
%             % nonlinear case is added in evaluate)
%             
%         end
    end
    
    methods(Access=protected)
        function f_eff = getf_eff(~, f, ~, ~)
            f_eff = f;
        end
        
        function B = getB_big(~, ~, C)
            null = zeros(size(C));
            B = [null, -eye(size(C)); null, C];
        end
    end
    
    methods(Access=private)
        function updatexJ(this, x)
            m = this.sys.Model;
            % Updates the x, J and R big versions (if needed)
            if isempty(this.x_big) || (norm(this.x_big - x) > 1e-8)
                [R_F, K_F] = this.weak_form(m, x);

                R = R_F(this.free);
                K = K_F(this.free, this.free);

                this.R_big = [zeros(size(R)); R];

                null = zeros(size(K));
                this.J_big = -1 * (this.B_big + [null, null; K, null]);
                this.x_big = x;
            end
        end
        
        function [R, K] = weak_form(this, m, u)
            % Nichtlineare Funktion (Schwache Form, R-f_s) und ihre Ableitung (nach den Freiheitsgraden) (tangentielle Steifigkeitsmatrix, K)

            data = m.data;
            
            % Steifigkeitsmatrix für u und T
            K = sparse(7 * data.num_knots, 7 * data.num_knots);
            % Residuumsvektor
            R = sparse(7 * data.num_knots, 1);

            % Tangentiale Steifigkeitsmatrix aufstellen und schwache Form mit der aktuellen Verschiebung auswerten
            el = this.sys.Model.Elements;
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