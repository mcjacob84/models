classdef DLTNonlinearCoreFun < dscomponents.ACoreFun
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
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties(Access=private)
        Kx;
        currentx = [];
        dKx;
    end
    
    methods
        function this = DLTNonlinearCoreFun(system)
            this = this@dscomponents.ACoreFun(system);
            
            free = this.System.Model.free;
            K0 = this.System.K0(free,free);
            this.fDim = size(K0,1);
            this.xDim = size(K0,2);
            this.JSparsityPattern = logical(K0);
        end
        
        function fx = evaluateCoreFun(this, x, t)%#ok
            % Daten nicht aktuell? (check inside)
            this.updateB(x);
            % Funktion f auswerten
            fx = -this.Kx;
        end
        
        function J = getStateJacobian(this, x, t, mu)%#ok
            % Daten nicht aktuell? (check inside)
            this.updateB(x);
            % Jacobi-Matrix ausgeben
            J = -this.dKx;
        end
        
        function pr = project(this, V, W)%#ok
            error('not implemented')
        end
    end
    
    methods(Access=private)
        function updateB(this, x)
            m = this.System.Model;
            % Updates the x, J and R big versions (if needed)
            if isempty(this.currentx) || (norm(this.currentx - x) > 1e-8)
                
                % Only pass the spatial part of x to the weak form!
                [Kx_full, dKx_full] = this.weak_form(x(1:length(m.free)));

                % K(x)
                this.Kx = Kx_full(m.free);
                
                % \partial K|\partial x
                % dKx = [0 0; dK 0]
                this.dKx = dKx_full(m.free, m.free);
                % Augment to represent zero velocity dof dependency
                this.dKx = [this.dKx zeros(size(this.dKx))];
                
                this.currentx = x;
            end
        end
        
        function [R, K] = weak_form(this, u)
            % Nichtlineare Funktion (Schwache Form, R-f_s) und ihre Ableitung (nach den Freiheitsgraden) (tangentielle Steifigkeitsmatrix, K)
            m = this.System.Model;
            data = m.data;
            % Full dimension
            fdim = 7*data.num_knots;
            % Steifigkeitsmatrix f�r u und T
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
                % sp�teren tests mal profilen. Idealerweise m�ssen die i,j
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
%             % Steifigkeitsmatrix f�r u und T
%             K = sparse(i,j,K,7 * data.num_knots, 7 * data.num_knots);
%             % Residuumsvektor
%             R = sparse(ir,ones(size(ir)),R, 7 * data.num_knots, 1);
        end
    end
end