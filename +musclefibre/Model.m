classdef Model < models.BaseFullModel
    % Model: Model for a muscle fibre model composed of motoneuron, spindle and
    % sarcomere array.
    %
    % The global time unit for this model is milliseconds [ms].
    %
    %
    % ToDos:
    % - Jacobimatrix-implementierung: evaluatePartialDerivatives. Dann könnten die eingebauten
    %  solver ggf schneller lösen, und der fehlerschätzer könnte implementiert werden
    % - Parameterbereich einschränken: nicht [0,1] sondern vielleicht erstmal [0,.2].
    % - Feststellen, welches der insgesamt zu betrachtende Simulationszeitraum sein sollte (abh.
    % vom aktuellen parameterbereich)
    % - Versuch, größere simulation für N>=400 oder so anzuwerfen
    % - Motoneuron untersuchen und parameterbereich auf sinnvolle größen einschränken! dazu eigenes
    % moto-model bauen und parameterbereich abtasten, dann parameterdomänen-geometrie einbauen.
    %
    % @author Daniel Wirtz @date 2012-11-22
    %
    % @new{0,8,dw,2015-09-15} Imported into +models package.
    %
    % @new{0,7,dw,2012-11-22} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.morepas.org/software/index.html
    % - \c Documentation http://www.morepas.org/software/kermor/index.html
    % - \c License @ref licensing
    
    properties
        % Seed that can be used by random number generator instances in order to enable result
        % reproduction.
        % @type int @default 1
        RandSeed = 1;
    end
    
    methods
        function this = Model(varargin)
            % Creates a new muscle fibre model
            %
            % Parameters:
            % N: number of sarcomeres @default 100
            
            i = inputParser;
            i.addParamValue('SarcoVersion',1);
            i.addParamValue('DynamicIC',true);
            i.addParamValue('SPM',false);
            i.addParamValue('OutputScaling',true);
            i.addParamValue('Spindle',false);
            i.addParamValue('N',100,@(n)n>1);
            i.parse(varargin{:});
            options = i.Results;
            
            this.Name = sprintf('Muscle fibre model (%d cells)',options.N);
            this.System = models.musclefibre.System(this, options);
            this.TrainingInputs = 1;
            this.SaveTag = sprintf('musclefibre_%dcells',options.N);
            this.Data = data.ModelData(this);
            if options.N > 1000
                this.Data.useFileTrajectoryData;
            end
            
            this.T = 400; % [ms]
            % DO NOT INCREASE! Peaks from Motoneuron are not resolved
            % visually correct if larger timesteps are used, leading to
            % confusion...
            this.dt = .1; % [ms]
            
            this.ODESolver = solvers.MLWrapper(@ode15s);
            
            this.Sampler.Domain = models.motoneuron.ParamDomain;
            
            this.DefaultMu = [.1; 3];
            this.DefaultInput = 1;
        end
        
        function pm = plot(this, t, y, pm)
            if nargin < 4
                pm = PlotManager(false,2,2);
                pm.LeaveOpen = true;
            end
            sys = this.System;
            off = sys.dm+sys.ds+(sys.MotoSarcoLinkIndex-1)*sys.dsa+1;
            h = pm.nextPlot('moto','Motoneuron V_s','time','value');
            plot(h,t,y(2,:));
            h = pm.nextPlot('sarco','Linked sarcomere: A_2','time','A_2');
            plot(h,t,y(off+52,:));
            h = pm.nextPlot('sarco','Linked sarcomere: V_s','time','V_s');
            plot(h,t,y(off,:));
            
            h = pm.nextPlot('sarco','All sarcomeres: V_s','time','value');
            Vsidx = sys.dm+sys.ds + 1:sys.dsa:size(y,1);
            plot(h,t,y(Vsidx,:));

            if nargin < 4
                pm.done;
            end    
        end
    end
    
    methods(Static, Access=protected)
        function this = loadobj(this)
            if ~isa(this, 'models.musclefibre.Model')
                sobj = this;
                this = models.musclefibre.Model(sobj.System.N);
                this.RandSeed = sobj.RandSeed;
                this = loadobj@models.BaseFullModel(this, sobj);
            else
                this = loadobj@models.BaseFullModel(this);
            end
        end
    end
    
    methods(Static)
        function res = test_MusclefibreModels
            res = true;
            for N = [2 5]
                for sv = [1 2]
                    for ic = logical([0 1])
                        try
                            fprintf('Testing N=%d,SV=%d,DynIC=%d\n',N,sv,ic);
                            m = models.musclefibre.Model('N',N,...
                                'SarcoVersion',sv,'DynamicIC',ic);
                            m.T = 40;
                            m.dt = .01;
                            m.simulate;

                            ms = models.musclefibre.Model('N',N,...
                                'SarcoVersion',sv,'DynamicIC',ic,...
                                'SPM',true);
                            ms.T = 40;
                            ms.dt = .01;
                            ms.simulate;
                        catch ME
                            display(ME)
                            res = false;
                            break;
                        end
                    end
                    if ~res
                        break;
                    end
                end
                if ~res
                    break;
                end
            end
        end
    end
end