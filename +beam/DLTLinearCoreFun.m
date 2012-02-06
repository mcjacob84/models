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
            % Set class for export
            this.CoeffClass = 'models.beam.dynlintimo.CoreFunCoeffs';
        end
        
        function initialize(this)
            initialize@models.beam.DLTBaseCoreFun(this);
            
            m = this.sys.Model;
            
            %% Fill inner AffLinCoreFun with matrices
            % Dämpfungsmodell 1: M a_t + (d1*M + d2*K) v_t + K u_t = f
            % d1 = 0.3 Dämpfungsfaktor vor Massenmatrix (Luftwiderstand)
            % d2 = .01 Dämpfungsfaktor vor Steifigkeitsmatrix (Materialdämpfung)
            % C = (mu(1)*M + mu(2)*K_0);
            K = this.sys.K0(m.free,m.free);
            s = length(m.free);
            null = sparse(s,s);
            % Clear affine parametric matrix
            this.AffParamMatrix.clear;
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