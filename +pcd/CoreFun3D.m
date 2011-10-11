classdef CoreFun3D < dscomponents.ACoreFun & ISimConstants
% The core nonlinear function of the PCD model in 3D.
%
% @author Daniel Wirtz @date 2011-10-05
%
% @new{0,5,dw,2011-10-05} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetObservable)
        % The number of seconds after which the neumann conditions decease.
        tDecaySecs = 2;
    end

    properties(Access=private)
        % The assoc. dynamical system
        % (need some values from that)
        sys;
        
        % The rect grid
        rbF;
        rbBa;
        rbL;
        rbR;
        rbT;
        rbBo;
        
        A;
        dim;
    end
    
    methods
        
        function this = CoreFun3D(dynsys)
            this.sys = dynsys;
            this.MultiArgumentEvaluations = false;
        end
        
        function copy = clone(this)
            copy = models.pcd.CoreFun3D(this.sys);
            
            % Call superclass method
            copy = clone@dscomponents.ACoreFun(this, copy);
            
            % copy reference!
            %copy.sys = this.sys; % already done in constructor
            copy.A = this.A;
            copy.dim = this.dim;
            copy.g = this.g;
        end
        
        function prepareConstants(this, mu, inputidx) %#ok
            % This fcn is called before each simulation is started.
            
            % Create diffusion matrix
            g = general.geometry.RectGrid3D(this.sys.dim1,this.sys.dim2,this.sys.dim3);
            this.A = general.MatUtils.laplacemat3D(this.sys.h, g);
            this.dim = g.Points;
            
            % Determine the points on each side that are affected of
            % boundary conditions (constant over each simulation as mu is
            % fixed then)
            h = this.sys.h;
            a = this.sys.Omega;
            xr = a(1,2)-a(1,1);
            yr = a(2,2)-a(2,1);
            zr = a(3,2)-a(3,1);
            %% Front & back (only x/y coords relevant)
            [i,j,k] = ind2sub(g.Dims,g.F);
            xd = abs((i-1)*h-.5*xr);
            yd = abs((j-1)*h-.5*yr);
            rb = xd < xr*mu(9)/2 & yd < yr*mu(9)/2;
            this.rbF = sub2ind(g.Dims, i(rb), j(rb), ones(size(i(rb))));
            rb = xd < xr*mu(10)/2 & yd < yr*mu(10)/2;
            this.rbBa = sub2ind(g.Dims, i(rb), j(rb), ones(size(i(rb)))*g.Dims(3));
            
            %% Left & Right (only x/z coords relevant)
            [i,j,k] = ind2sub(g.Dims,g.L);
            xd = abs((i-1)*h-.5*xr);
            zd = abs((k-1)*h-.5*zr);
            rb = xd < xr*mu(11)/2 & zd < zr*mu(11)/2;
            this.rbL = sub2ind(g.Dims, i(rb), ones(size(i(rb))), k(rb));
            rb = xd < xr*mu(12)/2 & zd < zr*mu(12)/2;
            this.rbR = sub2ind(g.Dims, i(rb), ones(size(i(rb)))*g.Dims(2), k(rb));
            
            %% Top & Bottom (only y/z coords relevant)
            [i,j,k] = ind2sub(g.Dims,g.T);
            yd = abs((j-1)*h-.5*yr);
            zd = abs((k-1)*h-.5*zr);
            rb = yd < yr*mu(13)/2 & zd < zr*mu(13)/2;
            this.rbT = sub2ind(g.Dims, ones(size(i(rb))), j(rb), k(rb));
            rb = yd < yr*mu(14)/2 & zd < zr*mu(14)/2;
            this.rbBo = sub2ind(g.Dims, ones(size(i(rb)))*g.Dims(1), j(rb), k(rb));
        end
        
        function fx = evaluateCoreFun(this, x, t, mu)
            % Allocate result vector
            fx = zeros(size(x));
            
            m = this.dim;
            n = this.sys.n;
            h = this.sys.h;
            
            if size(x,2) == 1
                % Extract single functions
                xa = x(1:m);
                xan = xa.^n;
                ya = x(m+1:2*m);
                xi = x(2*m+1:3*m);
                yi = x(3*m+1:end);

                % Compile boundary conditions
                rb = zeros(m,1);
                rb(this.rbF) = (xi(this.rbF)*mu(15)*Udecay(t))/h;
                rb(this.rbBa) = rb(this.rbBa) + (xi(this.rbBa)*mu(16)*Udecay(t))/h;
                rb(this.rbL) = rb(this.rbL) + (xi(this.rbL)*mu(17)*Udecay(t))/h;
                rb(this.rbR) = rb(this.rbR) + (xi(this.rbR)*mu(18)*Udecay(t))/h;
                rb(this.rbT) = rb(this.rbT) + (xi(this.rbT)*mu(19)*Udecay(t))/h;
                rb(this.rbBo) = rb(this.rbBo) + (xi(this.rbBo)*mu(20)*Udecay(t))/h;
%                 edges = [this.upper(1) this.upper(end)...
%                          this.lower(1) this.lower(end)];
%                 rb(edges) = .5*rb(edges);
                %rb(edges) = 0;

                fx(1:m) = mu(1)*xi.*ya - mu(3)*xa + this.A*xa + rb;
                fx(m+1:2*m) = mu(2)*yi.*xan - mu(4)*ya + this.sys.D2*this.A*ya;
                fx(2*m+1:3*m) = -mu(1)*xi.*ya - mu(5)*xi + mu(7) + this.sys.D3*this.A*xi - rb;
                fx(3*m+1:end) = -mu(2)*yi.*xan - mu(6)*yi + mu(8) + this.sys.D4*this.A*yi;
            else
                % Extract single functions
                xa = x(1:m,:);
                xan = xa.^n;
                ya = x(m+1:2*m,:);
                xi = x(2*m+1:3*m,:);
                yi = x(3*m+1:end,:);

                % Compile boundary conditions
                rb = zeros(m,size(x,2));
                rb(this.upper,:) = -(bsxfun(@mult,xi(this.upper,:),mu(9,:)))/h;
                rb(this.right,:) = rb(this.right,:) - (bsxfun(@mult,xi(this.right,:),mu(10,:)))/h;
                rb(this.lower,:) = rb(this.lower,:) - (bsxfun(@mult,xi(this.lower,:),mu(11,:)))/h;
                rb(this.left,:) = rb(this.left,:) - (bsxfun(@mult,xi(this.left,:),mu(12,:)))/h;
                edges = [this.upper(1) this.upper(end)...
                         this.lower(1) this.lower(end)];
                rb(edges,:) = .5*rb(edges,:);

                fx(1:m,:) = bsxfun(@mult,xi.*ya,mu(1,:)) - bsxfun(@mult,xa,mu(3,:)) + this.A*xa - rb;
                fx(m+1:2*m,:) = bsxfun(@mult,yi.*xan,mu(2,:)) - bsxfun(@mult,ya,mu(4,:)) + this.sys.D2*this.A*ya;
                fx(2*m+1:3*m,:) = -bsxfun(@mult,xi.*ya,mu(1,:)) - bsxfun(@mult,xi,mu(5,:)) + bsxfun(@mult,ones(size(xi)),mu(7,:)) + this.sys.D3*this.A*xi + rb;
                fx(3*m+1:end,:) = -bsxfun(@mult,yi.*xan,mu(2,:)) - bsxfun(@mult,yi,mu(6,:)) + bsxfun(@mult,ones(size(xi)),mu(8,:)) + this.sys.D4*this.A*yi;
            end
            
            function v = Udecay(t)
                if t < this.tDecaySecs
                    v = 1;
                else
                    v = max(0,2-.5*t);
                end
            end
            
            %         figure(1);
            %         plo(x,1);
            %         title('casp-8');
            %         plo(x,2);
            %         title('casp-3');
            %         plo(x,3);
            %         title('procasp-8');
            %         plo(x,4);
            %         title('procasp-3');
            %         pause;
        end
    end
    
%         function plo(fx,idx)
%         subplot(2,2,idx);
%         x = fx((idx-1)*m+1:idx*m);
%         surf(reshape(x,d1,d2));
%     end
end

