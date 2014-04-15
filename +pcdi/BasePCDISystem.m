classdef BasePCDISystem < models.BaseDynSystem
    %BasePCDISystem The base dynamical system class for the the Programmed
    % Cell Death Model by Markus Daub.
    %
    % See also: models.pcdi.BaseCoreFun
    %
    % @author Daniel Wirtz @date 2013-10-23
    %
    % @new{0,8,dw,2013-10-23} Added this class
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.morepas.org/software/index.html
    % - \c Documentation http://www.morepas.org/software/kermor/index.html
    % - \c License @ref licensing
    
    properties(SetObservable, Dependent)
        % Spatial stepwidth (in unscaled size units!)
        %
        % @propclass{critical} Determines the spatial resolution of the model.
        %
        % See also: L
        h; % is set in subclasses
        
        % The spatial width/area/region (in unscaled size units!)
        %
        % @propclass{data} The area for the model.
        %
        % See also: L
        Omega;
    end
    
    properties(Access=private)
        fh = [];
        fOmega = [];
    end
    
    properties(SetAccess=private)
        % Relative diffusion coefficients ([d2/d1, d3/d1, d4/d1])
        Diff;
        
        % scaled spatial stepwidth
        hs;
        
        % The system's dimensions
        Dims;
        
        % The reaction coefficients
        ReacCoeff;
        
        % The concentration labels
        Labels = {'Caspase-8','Caspase-3','Pro-Caspase-8','Pro-Caspase-3',...
            'IAP (iap)','BAR (bar)','Caspase-3+IAP (yb)','Caspase-8+BAR (xb)'};
        
        % The concentration image tags
        Tags = {'c8','c3','pc8','pc3','iap','bar','yb','xb'};
    end
    
    methods
        function this = BasePCDISystem(model)
            this = this@models.BaseDynSystem(model);
            
            % Scale diffusion coefficients
            m = this.Model;
            this.Diff = [m.d2/m.d1 m.d3/m.d1 m.d4/m.d1...
                         m.d5/m.d1 m.d6/m.d1 m.d7/m.d1 m.d8/m.d1];
            
            this.registerProps('h','Omega');
            
            this.ReacCoeff = [m.K1 m.K2 m.K3 m.K4 m.K5 m.K6 m.K7...
                              m.K8 m.K9 m.K10 m.K11 m.K12 m.K13...
                              m.Km3 m.Km8 m.Km9 m.Km10 m.Km11...
                              m.Km12]' * m.tau;
                          
            % Set state scaling
            this.StateScaling = this.Model.tc;
                          
            % Add parameters
            % Activation area
            p = this.addParam('area', [.01, 1], 10);
            p.Spacing = 'lin';
            
            % Activation rate
            p = this.addParam('rate', [0.0005, 1], 15);
            p.Spacing = 'log';
            
            % Activation time
            p = this.addParam('atime', [0, 1], 5);
            p.Spacing = 'lin';
            
            % Exponent
            p = this.addParam('exponent', [1, 4], 5);
            p.Spacing = 'lin';
            
            % Diffusion density
            p = this.addParam('diffusionstrength', [0, 1], 5);
            p.Spacing = 'lin';
        end
        
        function h = get.h(this)
            h = this.fh;
        end
        
        function set.h(this, value)
            if any(value >= this.fOmega(:,2))
                error('Cannot choose a step size h=%e value larger or equal to the geometry [%s].',...
                    value,Utils.implode(this.fOmega(:),', ','%g'));
            end
            this.fh = value;
            this.hs = value / this.Model.L;
            this.updateDims;
            
            m = this.Model;
            this.MaxTimestep = [];
            if ~isa(m.ODESolver,'solvers.IImplSolver') && ~isa(m.ODESolver,'solvers.SemiImplicitEuler')
                maxdt = .95*(value^2/max([m.d1 m.d2 m.d3 m.d4 m.d5 m.d6 m.d7 m.d8]));
                maxdtsc = .95*(this.hs^2/max([1 this.Diff]));
                
                if (maxdtsc - maxdt/this.Model.tau) /  maxdtsc > 1e-15
                    error('Inconsistent scaling. Please check.');
                end
                % Set max timestep value (scaled!)
                this.MaxTimestep = maxdtsc;
                
                if maxdt < m.dt
                    warning('models:pcd:CFL','CFL condition violated with h=%e and current dt=%e.\nSetting System.MaxTimestep=%e (scaled, effective value %e)\n', value, m.dt, maxdtsc, maxdt);
                end
            end
        end
        
        function v = get.Omega(this)
            v = this.fOmega;
        end
        
        function set.Omega(this, value)
            if size(value,2) ~= 2
                error('Omega needs to be a dim x 2 matrix with the spatial extend for each dimension in a row.');
            end
            this.fOmega = value;
            
            % Update the dimensions
            this.updateDims;
        end
    end
    
    methods(Access=private)
        
        function updateDims(this)
            if ~isempty(this.fh) && ~isempty(this.fOmega)
                nd = size(this.fOmega,1);
                this.Dims = zeros(1,nd);
                for d=1:nd
                    this.Dims(d) = length(this.fOmega(d,1):this.h:this.fOmega(d,2));
                end
                
                % Set new initial values and output in 1D-3D systems
                this.newSysDimension;
                
                % Call template method for custom actions in the core
                % functions (diffusion matrix comp)
                this.f.newSysDimension;
            end
        end
    end
    
    methods(Abstract, Access=protected)
        % Custom updates for new system dimension
        newSysDimension(this);
    end

end



