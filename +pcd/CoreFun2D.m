classdef CoreFun2D < dscomponents.ACoreFun
    % The core nonlinear function of the PCD model.
    %
    % @author Daniel Wirtz @date 16.03.2010
    %
    % @change{0,5,dw,2011-10-17} Removed the ISimConstants class and
    % unified the structure of the 1D-3D pcd models.
    
    properties(Access=private)
        % The assoc. dynamical system
        % (need some values from that)
        sys;
        
        idxmat;
        
        A;
        
        nodes;
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
        end
        
        function newSysDimension(this)
            % Create diffusion matrix
            [this.A, this.idxmat] = general.MatUtils.laplacemat(this.sys.h,...
                this.sys.Dims(1),this.sys.Dims(2));
            this.nodes = prod(this.sys.Dims);
        end
        
        function fx = evaluateCoreFun(this, x, t, mu)%#ok
            % Allocate result vector
            fx = zeros(size(x));
            
            m = this.nodes;
            s = this.sys;
            h = s.h;
            d = s.Dims;
            D = s.Diff;
            a = s.Omega;
            xr = a(1,2)-a(1,1);
            yr = a(2,2)-a(2,1);
            
            % Compile boundary conditions
            %to = this.sys.Model.tau*t;
            %ud = (to < 10)*1 + (to >= 10)*max(0,(2-to/10));
            ud = 1;
            
            if size(x,2) == 1
                % Extract single functions
                xa = x(1:m);
                xan = xa.^s.n;
                ya = x(m+1:2*m);
                xi = x(2*m+1:3*m);
                yi = x(3*m+1:end);
    
                rb = zeros(m,1);
                
                %% Top & Bottom
                xd = abs(((1:d(1))-1)*h-.5*xr); % y distances
                % bottom
                pos = xd < xr*mu(10)/2;
                idx = this.idxmat(pos,1); 
                rb(idx) = (xi(idx)*mu(14)*ud)/h;
                % top
                pos = xd < xr*mu(9)/2;
                idx = this.idxmat(pos,end);
                rb(idx) = rb(idx) + (xi(idx)*mu(13)*ud)/h;
                
                %% Left & Right
                yd = abs(((1:d(2))-1)*h-.5*yr);
                % right
                pos = yd < yr*mu(12)/2;
                idx = this.idxmat(end,pos);
                rb(idx) = rb(idx) + (xi(idx)*mu(16)*ud)/h;
                % left
                pos = yd < yr*mu(11)/2;
                idx = this.idxmat(1,pos);
                rb(idx) = rb(idx) + (xi(idx)*mu(15)*ud)/h;
                
                fx(1:m) = mu(1)*xi.*ya - mu(3)*xa + this.A*xa + rb;
                fx(m+1:2*m) = mu(2)*yi.*xan - mu(4)*ya + D(1)*this.A*ya;
                fx(2*m+1:3*m) = -mu(1)*xi.*ya - mu(5)*xi + mu(7) + D(2)*this.A*xi - rb;
                fx(3*m+1:end) = -mu(2)*yi.*xan - mu(6)*yi + mu(8) + D(3)*this.A*yi;
            else
                % Extract single functions
                xa = x(1:m,:);
                xan = xa.^s.n;
                ya = x(m+1:2*m,:);
                xi = x(2*m+1:3*m,:);
                yi = x(3*m+1:end,:);

                % Compile boundary conditions
                nd = size(x,2);
                rb = zeros(m,nd);
                
                %% Top & Bottom
                xd = repmat(abs(((1:d(1))-1)*h-.5*xr)',1,nd); % y distances
                % bottom
                pos = bsxfun(@lt,xd,xr*mu(10,:)/2);
                idx = this.idxmat(:,1);
                rb(idx,:) = pos .* (bsxfun(@mult,xi(idx,:),mu(14,:)*ud))/h;
                % top
                pos = bsxfun(@lt,xd,xr*mu(9,:)/2);
                idx = this.idxmat(:,end);
                rb(idx,:) = rb(idx,:) + pos .* (bsxfun(@mult,xi(idx,:),mu(13,:)*ud))/h;
                
                %% Left & Right
                yd = repmat(abs(((1:d(2))-1)*h-.5*yr)',1,nd);
                % right
                pos = bsxfun(@lt,yd,yr*mu(12,:)/2);
                idx = this.idxmat(end,:);
                rb(idx,:) = rb(idx,:) + pos .* (bsxfun(@mult,xi(idx,:),mu(16,:)*ud))/h;
                % left
                pos = bsxfun(@lt,yd,yr*mu(11,:)/2);
                idx = this.idxmat(1,:);
                rb(idx,:) = rb(idx,:) + pos .* (bsxfun(@mult,xi(idx,:),mu(15,:)*ud))/h;
                
                fx(1:m,:) = bsxfun(@mult,xi.*ya,mu(1,:)) - bsxfun(@mult,xa,mu(3,:)) + this.A*xa + rb;
                fx(m+1:2*m,:) = bsxfun(@mult,yi.*xan,mu(2,:)) - bsxfun(@mult,ya,mu(4,:)) + D(1)*this.A*ya;
                fx(2*m+1:3*m,:) = -bsxfun(@mult,xi.*ya,mu(1,:)) - bsxfun(@mult,xi,mu(5,:)) + bsxfun(@mult,ones(size(xi)),mu(7,:)) + D(2)*this.A*xi - rb;
                fx(3*m+1:end,:) = -bsxfun(@mult,yi.*xan,mu(2,:)) - bsxfun(@mult,yi,mu(6,:)) + bsxfun(@mult,ones(size(xi)),mu(8,:)) + D(3)*this.A*yi;
            end
        end
    end
end

