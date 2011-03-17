classdef BaseRBMatlabWrapper < models.BaseFullModel
    %BASERBMATLABMODEL Base class for all rbmatlab model wrappers
    %
    % Subclassing instances MUST explicitly call this classes constructor
    % in order to ensure that KerMor is connected to rbmatlab.
    %
    % @new{0,3,dw,2011-03-17} Added overloads for the
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
        RBMDataCont;
    end
    
    methods
        function this = BaseRBMatlabWrapper
            if ~KerMor.App.Hasrbmatlab
                error('rbmatlab is not registered with KerMor. Set KerMor.App.rbmatlabDirectory to fix.');
            end
        end
        
        function [t,y,sec,x] = simulate(this, mu, inputidx)
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
            [t,y,sec,x] = simulate@BaseModel(this, mu, inputidx);
        end
        
        function [t,x] = computeTrajectory(this, mu, inputidx)
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
            [t,x] = computeTrajectory@BaseModel(this, mu, inputidx);
        end
    end
    
end

