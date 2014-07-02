classdef Model < models.BaseFullModel
% MotoModel: Motoneuron model
%
% The global time unit for this model is milliseconds [ms].
%
% Relevante Frequenzbereiche: 6-50Hz.
%
% @todo place citations
%
% @todo Motoneuron untersuchen und parameterbereich auf sinnvolle gr��en einschr�nken! dazu eigenes
% moto-model bauen und parameterbereich abtasten, dann parameterdom�nen-geometrie einbauen.
%
% @author Daniel Wirtz @date 2012-11-22
%
% @new{0,7,dw,2013-07-05} Added this class.
%
    properties(Dependent)
        UseNoise;
    end
    
    methods
        function this = Model
            % Creates a new motoneuron model
            %
            this.Name = 'Motoneuron model';
            this.System = models.motoneuron.System(this);
            this.TrainingInputs = 1;
            this.SaveTag = 'motoneuron';
            this.Data = data.ModelData(this);
            
            this.T = 150; % [ms]
            this.dt = .1; % [ms]
            
            this.ODESolver = solvers.MLWrapper(@ode15s);
            
            this.DefaultMu = [.1; 3];
            this.DefaultInput = 1;
        end
        
        function varargout = plot(this, varargin)
            % plots some interesting states of the model
            %
            % See also: musclefibres.MuscleFibreSystem
            [varargout{1:nargout}] = this.System.plot(varargin{:});
        end
    end
    
    methods
        function value = get.UseNoise(this)
            value = ~this.System.noiseGen.DisableNoise;
        end
        
        function set.UseNoise(this, value)
            this.System.noiseGen.DisableNoise = ~value;
        end
    end
    
    methods(Static, Access=protected)
        function this = loadobj(this)
            if ~isa(this, 'models.motoneuron.Model')
                sobj = this;
                this = models.motoneuron.Model;
                this = loadobj@models.BaseFullModel(this, sobj);
            else
                this = loadobj@models.BaseFullModel(this);
            end
        end
    end
    
end