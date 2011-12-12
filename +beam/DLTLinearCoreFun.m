classdef DLTLinearCoreFun < models.beam.DLTBaseCoreFun & dscomponents.AffLinCoreFun
% DLTLinearCoreFun: 
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
    
    methods
        function this = DLTLinearCoreFun(system)
            this = this@models.beam.DLTBaseCoreFun(system);
            this.initialize;
        end
        
        function initialize(this)
            initialize@models.beam.DLTBaseCoreFun(this);
            
            m = this.sys.Model;
            
            %% Add constant dirichlet forces
            % Reconstruct fake u vector which has entries only at
            % dirichlet points. Used for computation of constant forces due
            % to dirichlet values.
            u = zeros(7*m.data.num_knots,1);
            u(m.dir_u) = m.u_dir;
            f_dir = -this.K0*u; 
            d = length(m.free);
            this.f_big(d+1:end) = this.f_big(d+1:end) + f_dir(m.free);
            
            %% Assign offset part of linear core function
            this.b = this.f_big;
            
            %% Fill inner AffLinCoreFun with matrices
            % Dämpfungsmodell 1: M a_t + (d1*M + d2*K) v_t + K u_t = f
            % d1 = 0.3 Dämpfungsfaktor vor Massenmatrix (Luftwiderstand)
            % d2 = .01 Dämpfungsfaktor vor Steifigkeitsmatrix (Materialdämpfung)
            % C = (mu(1)*M + mu(2)*K);
            K = this.K0(m.free,m.free);
            s = length(m.free);
            null = sparse(s,s);
            % Add constant term
            this.addMatrix('1',[null -speye(s);...
                                K null]);
                            % Remember: this.M is a ConstantMassMatrix
                            % class which has the M property as matrix!
            this.addMatrix('mu(1)',[null null;...
                                    null this.sys.M_small]);
            this.addMatrix('mu(2)',[null null;...
                                    null K]);
        end
    end
end