classdef System < models.BaseFirstOrderSystem
    % System: The global dynamical system used within the MotoModel
    %
    %
    % @author Timm Strecker @date 2014-03-24
    %
    % @new{0,7,ts,2014-03-24} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    
    properties(Transient)
        % debug output
        d = [];
    end
    
    
    properties(Dependent)
        % the conductivity tensor in intracellular space
        %
        % @type matrix<double>
        sigmai;
        
        % the conductivity tensor in extracellular space
        %
        % @type matrix<double>
        sigmae;
        
        % the conductivity tensor in the body (skin/fat)
        %
        % @type matrix<double>
        sigmao;
    end
    
    properties(Access=private)
        % the conductivity tensor in intracellular space
        %
        % @type matrix<double>
        sigmai_ = diag([8.93 0.893 0.893]); %eye(3);
        
        % the conductivity tensor in extracellular space
        %
        % @type matrix<double>
        sigmae_ = diag([6.7 6.7 6.7]);  %eye(3);
        
        % the conductivity tensor in the body (skin/fat)
        %
        % @type matrix<double>
        sigmao_ = diag([.4 .4 .4]);  %eye(3);
    end
    
    properties(SetAccess = private)
        % Size of the muscle. vector with 4 entries representing length,
        % width, muscle height and height of the fat/skin layer. in cm
        %
        % @default [1;1;1;0]
        musclegeometry;
        
        % Discretization stepwidth in cm.
        %
        % @default [0.1; 0.1; 0.1]
        h;
        
        % Dimension of the discretization.
        %
        % @default [10;10;10]
        dim;
        
        % Dimension of the input u
        %
        % @default 1000
        udim;
        
        % number of nodes in z direction lying within the muscle domain
        %
        % @default 10
        zdim_m;
    end
    
    methods
        function this = System(model, musclegeometry, dim)
            % Call superclass constructor
            this = this@models.BaseFirstOrderSystem(model);
            
            this.updateDimensions(musclegeometry, dim);
            
            % Parameters mu
            this.addParam('mean current', 4, 'Range', [0 9], 'Desired', 100); 
            
            this.f = models.emg.RHS(this);
            this.assembleA;
        end
        
        function pm = plot(this, t, y, pm)
            if nargin < 4
                pm = PlotManager;
                pm.LeaveOpen = true;
            end
            xaxis = 0:this.h(1):this.musclegeometry(1);
            yaxis = 0:this.h(2):this.musclegeometry(2);
            ax = pm.nextPlot('emg','','width [cm]','length [cm]');
            zlabel(ax,'\phi_i [mV]');
            pm.done;
            view(53,28);
            surfacepos = reshape(1:size(y,1),this.dim(1),this.dim(2),[]);
            surfacepos = surfacepos(:,:,end);
            hlp = y(surfacepos(:),:)*10;
            %axis([0 this.musclegeometry(2) 0 this.musclegeometry(1) min(hlp(:)) max(hlp(:))]);
            maxval = ceil(max(max(y(surfacepos,:))));
            minval = floor(min(min(y(surfacepos,:))));
            axis([0 this.musclegeometry(2) 0 this.musclegeometry(1) minval maxval]);
            caxis([minval maxval]);
            % daspect([1 1 1]);
            daspect([1 1 (maxval-minval)]);
            hold(ax,'on');
            for idx = 1:length(t)
                cla(ax);
                phi = reshape(y(:,idx),this.dim(1),this.dim(2),[]);
                % surf(ax, yaxis, xaxis, 10*phi(:,:,end),'FaceColor','interp','EdgeColor','interp');  % pcolor?
                surf(ax, yaxis, xaxis, phi(:,:,end),'FaceColor','interp','EdgeColor','interp');  % pcolor?
                title(ax,sprintf('Surface EMG at time %gms',t(idx)));
                colorbar;
                pause(.2);
                if ~ishandle(ax)
                    break;
                end
            end
            if nargin < 4
                pm.done;
            end
        end
        
        function assembleA(this)
            ns = this.NumStateDofs;
            musclegrid = general.geometry.RectGrid3D(this.dim(1),this.dim(2),this.dim(3));
            [ie, je, se] = find(MatUtils.generalizedLaplacemat3D(this.h,musclegrid,this.sigmai_+this.sigmae_));
            io = []; jo = []; so = [];
            ib = []; jb = []; sb = [];  % not used if A is constructed this way
            
            % this variant works best:
            if this.musclegeometry(4) ~= 0  % there is a skin layer
                % strategy: create grids for muscle and skin that are too
                % high. create corresponding generalized Laplacians. remove
                % boundary elements between muscle and skin and concatenate Laplacians.

                idxmax = musclegrid.IndexMatrix(end,end,this.zdim_m);
                idx = (ie<=idxmax);
                %ie = ie(1:mmax); je = je(1:mmax); se = se(1:mmax);
                ie = ie(idx); je = je(idx); se = se(idx);
                skingrid = general.geometry.RectGrid3D(this.dim(1),this.dim(2),this.dim(3)-this.zdim_m+1);
                [io, jo, so] = find(MatUtils.generalizedLaplacemat3D(this.h,skingrid,this.sigmao_));
                idxmin = skingrid.IndexMatrix(1,1,2);
                idx = (io >= idxmin);
                %io = io(smin:end); jo = jo(smin:end); so = so(smin:end);
                io = io(idx); jo = jo(idx); so = so(idx);
                io = io + this.udim - idxmin + 1; jo = jo + this.udim - idxmin + 1;
                % ToDo: consider d/dx sigma (because of product rule) - but not
                % differentiable...
            end
         
            % integral over domain
            ii = (ns+1) * ones(1,ns);
            ji = 1:ns;
            si = prod(this.h)*ones(1,ns);
            
            A = sparse([ie' ib' io' ii], [je' jb' jo' ji], [se' -sb' so' si], ns+1,ns);
            this.A = dscomponents.LinearCoreFun(A);
        end        

    end
    
    methods(Access=protected)
        function updateDimensions(this, musclegeometry, dim)
            this.dim = dim;
            this.zdim_m = round(musclegeometry(3)/(musclegeometry(3)+musclegeometry(4))*dim(3));
            this.musclegeometry = musclegeometry;
            this.udim = dim(1)*dim(2)*this.zdim_m;
            this.h = [musclegeometry(1:2);musclegeometry(3)+musclegeometry(4)]./(dim-1);
            
            this.NumStateDofs = prod(dim);
            updateDimensions@models.BaseFirstOrderSystem(this);
        end
    end
    
    methods %% setters
        
        function set.sigmae(this, value)
            this.sigmae_ = value;
            this.assembleA;
        end
        
        function value = get.sigmae(this)
            value = this.sigmae_;
        end
        
        function set.sigmai(this, value)
            this.sigmai_ = value;
            this.assembleA;
            this.f.assembleB;
        end
        
        function value = get.sigmai(this)
            value = this.sigmai_;
        end
        
        function set.sigmao(this, value)
            this.sigmao_ = value;
            this.assembleA;
        end
        
        function value = get.sigmao(this)
            value = this.sigmao_;
        end
        
    end
end

