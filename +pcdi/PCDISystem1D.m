classdef PCDISystem1D < models.pcdi.BasePCDISystem
    %PCDISystem1D Inhibited PCD system 1D variant
    %
    % @author Daniel Wirtz @date 2013-10-23
    %
    % @new{0,8,dw,2013-10-23} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    
    methods
        function this = PCDISystem1D(model)
            this = this@models.pcdi.BasePCDISystem(model);

            % Set core function
            this.f = models.pcdi.CoreFun1D(this);
            
            % Spatial resolution (in real sizes)
            this.Omega = [0 1] * this.Model.L;
            this.h = this.Model.L/24;
           
            % Remove params set in BasePCDSystem (simpler for 1D)
            this.Params = data.ModelParam.empty;
            
            % Add input param (is getting inserted after the BasePCDSystem
            % constructor params, so number 9!)
            this.addParam('U', [1e-5, 1e-2], 50);
        end

        function plot(this, model, t, y)
            % Performs a plot for this model's results.
            %
            
            if length(t) > 700
                idx = round(linspace(1,length(t),700));
                t = t(idx);
                y = y(:,idx);
            end
            states = {'alive','unstable','dead'};
            
            figure;
            X = t;
            Y = this.Omega(1):this.h:this.Omega(2);
            m = model.System.Dims(1);
            doplot(y(1:m,:),'Caspase-8 (x_a)',1);
            doplot(y(m+1:2*m,:),'Caspase-3 (y_a)',2);
            doplot(y(2*m+1:3*m,:),'Pro-Caspase-8 (x_i)',3);
            doplot(y(3*m+1:end,:),'Pro-Caspase-3 (y_i)',4);
            
            function doplot(y,thetitle,pnr)
                subplot(2,2,pnr);
                mesh(X,Y,y);
                %surf(X,Y,y,'EdgeColor','none');
                xlabel('Time [s]');
                ylabel(sprintf('%.2e to %.2e: Cell core to hull [m]',this.Omega(1),this.Omega(2)));
                grid off;
                mi = min(y(:));
                Ma = max(y(:));
                if abs((mi-Ma) / mi) < 1e-14
                    mi = .999*mi; Ma=1.001*Ma;
                end
                axis([0 max(t) this.Omega mi Ma]);
                di = abs(this.SteadyStates(:,pnr)-y(end));
                reldi = di ./ (this.SteadyStates(:,pnr)+eps);
                reldistr = Utils.implode(reldi,', ','%2.3e');
                if any(reldi > .1) || any(reldi < 10)
                    [~, id] = min(di);
                    title(sprintf('Model "%s", %s concentrations\nCell state at T=%d: %s\n%s', model.Name, thetitle,...
                    max(t),states{id},reldistr));
                else
                    title(sprintf('Model "%s", %s concentrations\n%s', model.Name, thetitle,reldistr));
                end
            end
        end
    end
    
    methods(Access=protected)
        function newSysDimension(this)
            % Assign fitting initial value
            m = prod(this.Dims);
            
            x0 = zeros(4*m,1);
            x0(1:2*m) = 1e-16;
            x0(2*m+1:end) = 1e-9;
            this.x0 = dscomponents.ConstInitialValue(x0);
            
            e = ones(m,1);
            A = spdiags([e -2*e e],-1:1,m,m)/this.hs^2;
            A(1) = -1/this.hs^2;
            A(end) = -1/this.hs^2;
            A = blkdiag(A,this.Diff(1)*A,this.Diff(2)*A,this.Diff(3)*A);
            this.A = dscomponents.LinearCoreFun(A);
            
            % Extracts the caspase-3 concentrations from the result
%             C = zeros(m,4*m);
%             C(:,m+1:2*m) = diag(ones(m,1));
            %this.C = dscomponents.LinearInputConv(sparse(C));
        end
    end
end

