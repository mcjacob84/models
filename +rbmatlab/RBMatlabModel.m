classdef RBMatlabModel < models.BaseFullModel
    %RBMatlabModel: Base class for all rbmatlab models in KerMor
    %
    % Subclassing instances MUST explicitly call this classes constructor
    % in order to ensure that KerMor is connected to rbmatlab.
    %
    % @change{0,3,sa,2011-05-11} Implemented setters for the properties
    %
    % @new{0,2,dw,2011-03-17} Added overloads for the
    % models.BaseModel.simulate and models.BaseModel.computeTrajectory
    % methods that check for a valid connection to rbmatlab before
    % computations are made.
    
    properties
        % The adopted RB matlab model.
        %
        % Simply wraps the variable 'model'.
        RBMModel;
        
        % The RBMDataContainer for the model_data struct
        %
        % A handle class containing the usual 'model_data' struct as
        % single property RBMData.
        %
        % @type models.rbmatlab.RBMDataContainer
        RBMDataCont;
    end
    
    methods
        function this = RBMatlabModel
            if ~KerMor.App.Hasrbmatlab
                error('rbmatlab is not registered with KerMor. Set KerMor.App.rbmatlabDirectory to fix.');
            end
            this = this@models.BaseFullModel;
        end
        
        function set.RBMModel(this, value)
            if ~isstruct(value)
                error('Value must be a valid double matrix');
            end
            this.RBMModel = value;
        end
        
        function set.RBMDataCont(this, value)
            this.checkType(value, 'models.rbmatlab.RBMDataContainer');
            this.RBMDataCont = value;
        end
        
        function [t,y,sec,x] = simulate(this, varargin)
            % Simulates the RBMatlab model for given parameter and
            % inputindex.
            % 
            % Overloads the base method in models.BaseModel.
            %
            % Inherited documentation:
            % @copydoc models::BaseModel::simulate()
            if ~KerMor.App.Hasrbmatlab
                error('rbmatlab is not registered with KerMor. Set KerMor.App.rbmatlabDirectory to fix.');
            end
            [t,y,sec,x] = simulate@models.BaseModel(this, varargin{:});
        end
        
        function [t,x] = computeTrajectory(this, varargin)
            % Computes a trajectory of the RBMatlab model for given
            % parameter and inputindex.
            % 
            % Overloads the base method in models.BaseModel.
            %
            % Inherited documentation:
            % @copydoc models::BaseModel::computeTrajectory()
            if ~KerMor.App.Hasrbmatlab
                error('rbmatlab is not registered with KerMor. Set KerMor.App.rbmatlabDirectory to fix.');
            end
            [t,x] = computeTrajectory@models.BaseModel(this, varargin{:});
        end
    end
    
end

