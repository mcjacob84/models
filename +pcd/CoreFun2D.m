classdef CoreFun2D < dscomponents.ACoreFun
% The core nonlinear function of the PCD model.
%
% @author Daniel Wirtz @date 16.03.2010
%
% @change{0,5,dw,2011-11-02} Augmenting the mu parameters by the base system's
% models.pcd.BasePCDSystem.ReacCoeff vector. This removes the reaction coefficients from
% the system as true parameters but allows to quickly revert the process if needed.
%
% @change{0,5,dw,2011-10-17} Removed the ISimConstants class and
% unified the structure of the 1D-3D pcd models.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(Access=private)
        % The assoc. dynamical system
        % (need some values from that)
        sys;
        
        idxmat;
        
        A;
        
        nodes;
        
        % Contains helper values that stay constant for each geometry
        %
        % See newSysDimension for details
        hlp = struct;
    end
    
    methods
        
        function copy = clone(this)
            % sys already copied in constructor (see below)
            copy = models.pcd.CoreFun2D(this.sys);
            
            % Call superclass method
            copy = clone@dscomponents.ACoreFun(this, copy);
            
            copy.A = this.A;
            copy.idxmat = this.idxmat;
            copy.nodes = this.nodes;
        end
        
        function this = CoreFun2D(dynsys)
            this.sys = dynsys;
            this.MultiArgumentEvaluations = true;
            
            this.hlp.Diff = this.sys.Diff;
            this.hlp.n = this.sys.n;
        end
        
        function newSysDimension(this)
            % Create diffusion matrix
            s = this.sys;
            [this.A, this.idxmat] = general.MatUtils.laplacemat(s.hs,...
                s.Dims(1),s.Dims(2));
            this.nodes = prod(s.Dims);
            this.hlp.d1 = s.Dims(1);
            this.hlp.d2 = s.Dims(2);
            
            % Can use the unscaled h and original Omega for computation of
            % ranges & distances from middle
            this.hlp.hs = s.hs;
            a = s.Omega;
            this.hlp.xr = a(1,2)-a(1,1);
            this.hlp.yr = a(2,2)-a(2,1);
            this.hlp.xd = abs(((1:this.hlp.d1)-1)*this.sys.h-.5*this.hlp.xr); % x distances
            this.hlp.yd = abs(((1:this.hlp.d2)-1)*this.sys.h-.5*this.hlp.yr); % y distances
        end
        
        function fx = evaluateCoreFun(this, x, t, mu)%#ok
            % Allocate result vector
            fx = zeros(size(x));
            
            m = this.nodes;
            
            % Compile boundary conditions
            %to = this.sys.Model.tau*t;
            %ud = (to < 10)*1 + (to >= 10)*max(0,(2-to/10));
            ud = 1;
            
            if size(x,2) == 1
                % Uncomment if reaction coeffs become real params again
                %mu = [s.ReacCoeff; mu]';
                mu = [this.sys.ReacCoeff; mu([1 1 1 1 2 2 2 2])];
                
                % Extract single functions
                xa = x(1:m);
                xan = xa.^this.hlp.n;
                ya = x(m+1:2*m);
                xi = x(2*m+1:3*m);
                yi = x(3*m+1:end);
    
                rb = zeros(m,1);
                
                %% Top & Bottom
                % bottom
                pos = this.hlp.xd <= this.hlp.xr*mu(10)/2;
                idx = this.idxmat(pos,1); 
                rb(idx) = (xi(idx)*mu(14)*ud);
                % top
                pos = this.hlp.xd <= this.hlp.xr*mu(9)/2;
                idx = this.idxmat(pos,end);
                rb(idx) = rb(idx) + (xi(idx)*mu(13)*ud);
                
                %% Left & Right
                % right
                pos = this.hlp.yd <= this.hlp.yr*mu(12)/2;
                idx = this.idxmat(end,pos);
                rb(idx) = rb(idx) + (xi(idx)*mu(16)*ud);
                % left
                pos = this.hlp.yd <= this.hlp.yr*mu(11)/2;
                idx = this.idxmat(1,pos);
                rb(idx) = rb(idx) + (xi(idx)*mu(15)*ud);
                
                fx(1:m) = mu(1)*xi.*ya - mu(3)*xa + this.A*xa + rb/this.hlp.hs;
                fx(m+1:2*m) = mu(2)*yi.*xan - mu(4)*ya + this.hlp.Diff(1)*this.A*ya;
                fx(2*m+1:3*m) = -mu(1)*xi.*ya - mu(5)*xi + mu(7) + this.hlp.Diff(2)*this.A*xi - rb/this.hlp.hs;
                fx(3*m+1:end) = -mu(2)*yi.*xan - mu(6)*yi + mu(8) + this.hlp.Diff(3)*this.A*yi;
            else
                % Uncomment if reaction coeffs become real params again
                %mu = [s.ReacCoeff; mu]';
                mu = [repmat(this.sys.ReacCoeff,1,size(mu,2)); mu([1 1 1 1 2 2 2 2],:)];
                
                % Extract single functions
                xa = x(1:m,:);
                xan = xa.^this.hlp.n;
                ya = x(m+1:2*m,:);
                xi = x(2*m+1:3*m,:);
                yi = x(3*m+1:end,:);

                % Compile boundary conditions
                nd = size(x,2);
                rb = zeros(m,nd);
                
                %% Top & Bottom
                xd = repmat(this.hlp.xd',1,nd); % y distances
                % bottom
                pos = bsxfun(@lt,xd,this.hlp.xr*mu(10,:)/2);
                idx = this.idxmat(:,1);
                rb(idx,:) = pos .* (bsxfun(@mult,xi(idx,:),mu(14,:)*ud));
                % top
                pos = bsxfun(@lt,xd,this.hlp.xr*mu(9,:)/2);
                idx = this.idxmat(:,end);
                rb(idx,:) = rb(idx,:) + pos .* (bsxfun(@mult,xi(idx,:),mu(13,:)*ud));
                
                %% Left & Right
                yd = repmat(this.hlp.yd',1,nd);
                % right
                pos = bsxfun(@lt,yd,this.hlp.yr*mu(12,:)/2);
                idx = this.idxmat(end,:);
                rb(idx,:) = rb(idx,:) + pos .* (bsxfun(@mult,xi(idx,:),mu(16,:)*ud));
                % left
                pos = bsxfun(@lt,yd,this.hlp.yr*mu(11,:)/2);
                idx = this.idxmat(1,:);
                rb(idx,:) = rb(idx,:) + pos .* (bsxfun(@mult,xi(idx,:),mu(15,:)*ud));
                
                fx(1:m,:) = bsxfun(@mult,xi.*ya,mu(1,:)) - bsxfun(@mult,xa,mu(3,:)) + this.A*xa + rb/this.hlp.hs;
                fx(m+1:2*m,:) = bsxfun(@mult,yi.*xan,mu(2,:)) - bsxfun(@mult,ya,mu(4,:)) + this.hlp.Diff(1)*this.A*ya;
                fx(2*m+1:3*m,:) = -bsxfun(@mult,xi.*ya,mu(1,:)) - bsxfun(@mult,xi,mu(5,:)) + bsxfun(@mult,ones(size(xi)),mu(7,:)) + this.hlp.Diff(2)*this.A*xi - rb/this.hlp.hs;
                fx(3*m+1:end,:) = -bsxfun(@mult,yi.*xan,mu(2,:)) - bsxfun(@mult,yi,mu(6,:)) + bsxfun(@mult,ones(size(xi)),mu(8,:)) + this.hlp.Diff(3)*this.A*yi;
            end
        end
    end
end

