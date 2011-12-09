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
            for i = 1:data.num_elem_RO
                index_0_glob = 7*data.knot_index(m.RO(i).PointsIdx(1))-6 : 7*data.knot_index(m.RO(i).PointsIdx(1))-1;
                index_L_glob = 7*data.knot_index(m.RO(i).PointsIdx(2))-6 : 7*data.knot_index(m.RO(i).PointsIdx(2))-1;
                index_glob = [index_0_glob index_L_glob];

                u_e = u(index_glob);

                [K_lok, R_lok] = m.RO(i).getLocalTangentials(u_e);

                K(index_glob, index_glob) = K(index_glob, index_glob) + K_lok;
                R(index_glob) = R(index_glob) + R_lok;
            end

            for i = 1:data.num_elem_KR
                index_0_glob = 7*data.knot_index(m.KR(i).PointsIdx(1))-6 : 7*data.knot_index(m.KR(i).PointsIdx(1))-1;
                index_L_glob = 7*data.knot_index(m.KR(i).PointsIdx(2))-6 : 7*data.knot_index(m.KR(i).PointsIdx(2))-1;
                index_glob = [index_0_glob index_L_glob];

                u_e = u(index_glob);

                [K_lok, R_lok] = m.KR(i).getLocalTangentials(u_e);

                K(index_glob, index_glob) = K(index_glob, index_glob) + K_lok;
                R(index_glob) = R(index_glob) + R_lok;
            end

            for i = 1:data.num_elem_FH
                index_0_glob = 7*data.knot_index(m.FH(i).PointsIdx(1))-6 : 7*data.knot_index(m.FH(i).PointsIdx(1))-4;
                index_L_glob = 7*data.knot_index(m.FH(i).PointsIdx(2))-6 : 7*data.knot_index(m.FH(i).PointsIdx(2))-4;
                index_glob = [index_0_glob index_L_glob];

                u_e = u(index_glob);

                [K_lok, R_lok] = m.FH(i).getLocalTangentials(u_e);

                K(index_glob, index_glob) = K(index_glob, index_glob) + K_lok;
                R(index_glob) = R(index_glob) + R_lok;
            end
        end
    end
end