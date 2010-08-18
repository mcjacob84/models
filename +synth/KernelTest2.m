classdef KernelTest2 < models.BaseFullModel & models.BaseDynSystem & dscomponents.ACoreFun
    % Kernel core function test model 1
    % 
    % This class implements both the model and dynamical system!
    %
    
    properties
        tk;
        pk;
        sk;
        
        % The system's dimension
        dim = 500;
        
        % SV's
        sv = 10;
    end
    
    properties(Access=private)
        xi;
        mui;
        ti;
        ai;
    end
    
    methods
        function this = KernelTest1(dims)
            
            if nargin == 0
                dims = 500;
            end
            
            %% Model settings
            this.Verbose = 0;
            this.Name = 'Synthetic kernel based test model 1';
            
            this.T = 5;
            this.dt = .1;
            
            this.Sampler = sampling.GridSampler;
            
            %a = approx.CompWiseInt;
            %a.TimeKernel = kernels.LinearKernel;
            %a.SystemKernel = kernels.GaussKernel(.5);
            %a.ParamKernel = kernels.GaussKernel(.5);
            %this.Approx = a;
            
            s = spacereduction.PODReducer;
            s.Mode = 'abs';
            s.Value = 1;
            this.SpaceReducer = s;
            
            this.ODESolver = @ode45;
            
            %% System settings
            this.f = this;
            
             % Add params
            %this.addParam('x0 multiplier', [-1 1], 4);
            %this.addParam('Random param 1', [-1 1], 5);
            %this.addParam('Random param 2', [-3 3], 10);
            
            this.Inputs{1} = @(t)t/2;
            
            % Assign same kernels
            s.sk = a.SystemKernel;
            this.System = s;
        end
        
        function fx = evaluate(this, x, t, mu)
%             fx = this.ai * (this.tk.evaluate(this.ti,t)...
%                 .* this.pk.evaluate(this.mui,mu)...
%                 .* this.sk.evaluate(this.xi, x));
            % only x is used!
            fx = this.ai * this.sk.evaluate(this.xi, x);
            fx = repmat(fx,this.dim,1);
        end
        
        function c = getGlobalLipschitz(this, t, mu)
            c = norm(this.ai) * this.sk.getGlobalLipschitz;
        end
        
        function projected = project(this, V)
            % Implements ACoreFun(::IProjectable).project
            % no projection for this model as there are only two
            % dimensions.
            projected = this;
        end
        
        function showBaseFun(this)
            % Debug method; displays the core function for each parameter
            % sample.
            %for muidx = 1:size(this.mui,2)
                figure(1);
                
                
                %t = 0:this.dt/2:this.T;
                %numt = length(t);
                %x = linspace(-10,10,this.sv*5);
                %allx = x;
                %numx = size(x,2);
                
                %allx = repmat(x,1,numt);
                %repidx = repmat((1:numt)',1,numx);
                %tidx = reshape(repidx',1,[]);
                %allt = t(tidx);
                
                % Fix one parameter
                %mu = repmat(this.mui(:,muidx),1,size(allx,2));
                
%                 fx = this.ai * (this.tk.evaluate(this.ti,allt)...
%                     .* this.pk.evaluate(this.mui, mu)...
%                     .* this.sk.evaluate(this.xi(1,:),allx));
                x = repmat(linspace(-10,10,this.sv*5),this.dim,1);
                fx = this.evaluate(x);
%                 fx1 = this.ai * this.sk.evaluate(this.xi(1,:),allx);
%                 fx2 = this.ai * this.sk.evaluate(this.xi(1:2,:),[allx; allx]);
%                 fx3 = this.ai * this.sk.evaluate(this.xi(1:3,:),[allx; allx; allx]);
                plot(x(1,:),fx(1,:),'r');
                %plot(x,fx1,'r',x,fx2,'b',x,fx3,'g');
                
%                 fx = reshape(fx,numx,numt);
%                 [X,Y] = meshgrid(t,x);
%                 surfl(X,Y,fx);
%                 axis tight;
%                 xlabel('t');
%                 ylabel('x');
%                 s = sprintf('Plot for parameter mu=%f',mu);
%                 title(s);
%                 pause;
            %end
        end
        
    end
    
end

