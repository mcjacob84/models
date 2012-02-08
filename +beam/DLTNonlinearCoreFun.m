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
        Kx_big;
        currentx = [];
        B_big;
        B_big_const;
        dKx_big;
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
            fx = -this.B_big*x - this.Kx_big;
        end
        
        function J = getStateJacobian(this, x, t, mu)%#ok
            % Daten nicht aktuell? (check inside)
            this.updateB(x, mu);
            % Jacobi-Matrix ausgeben
            J = -this.B_big - this.dKx_big;
        end
        
        function pr = project(this, V, W)%#ok
            error('not implemented')
        end
    end
    
    methods(Access=private)
        function updateB(this, x, mu)
            m = this.sys.Model;
            % Updates the x, J and R big versions (if needed)
            if isempty(this.currentx) || (norm(this.currentx - x) > 1e-8)
                s = length(m.free);
                
                % Only pass the spatial part of x to the weak form!
                [Kx, dKx] = this.weak_form(x(1:length(m.free)));

                % K(x)
                Kx = Kx(m.free);
                this.Kx_big = [0*Kx; Kx];
                
                % \partial K|\partial x
                % dKx = [0 0; dK 0]
                this.dKx_big = sparse(2*s,2*s);
                this.dKx_big(s+1:end,1:s) = dKx(m.free, m.free);
                
                % Assemble damping matrix (dep. on parameter)
                K0 = this.sys.K0(m.free,m.free);
                C_big = sparse(2*s,2*s);
                C_big(s+1:end,s+1:end) = mu(1)*K0 + mu(2)*this.sys.M_small;
                % B = [0 -I; 0 C]
                this.B_big = this.B_big_const + C_big;
                
                this.currentx = x;
            end
        end
        
        function [R, K] = weak_form(this, u)
            % Nichtlineare Funktion (Schwache Form, R-f_s) und ihre Ableitung (nach den Freiheitsgraden) (tangentielle Steifigkeitsmatrix, K)
            m = this.sys.Model;
            data = m.data;
            % Full dimension
            fdim = 7*data.num_knots;
            % Steifigkeitsmatrix für u und T
            K = sparse(fdim, fdim);
            % Residuumsvektor
            R = sparse(fdim, 1);

            % Tangentiale Steifigkeitsmatrix aufstellen und schwache Form mit der aktuellen Verschiebung auswerten
            el = m.Elements;
%             i = []; j = []; K = []; R = []; ir = [];
            u_full = zeros(fdim,1);
            % load dirichlet values
            u_full(m.dir_u) = m.u_dir;
            % set current values, only use x
            u_full(m.free) = u;
            for k=1:length(el)
                index_glob = el{k}.getGlobalIndices;
                [K_lok, R_lok] = el{k}.getLocalTangentials(u_full(index_glob));
                
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