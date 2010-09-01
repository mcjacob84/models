classdef BaseRBMatlabWrapper < models.BaseFullModel
    %BASERBMATLABMODEL Base class for all rbmatlab model wrappers
    %
    % Subclassing instances MUST explicitly call this classes constructor
    % in order to gain access to rbmatlab's file paths etc.
    
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
            % Environmental setup
            rbmatlab = '/afs/.mathe/project/agh/home/dwirtz/rbmatlab/';
            if isempty(getenv('KERMORTEMP'))
                error('No KERMORTEMP environment path set. Aborting. (Forgot to run startup_kermor?)');
            end
            setenv('RBMATLABTEMP', getenv('KERMORTEMP'));
            setenv('RBMATLABHOME',rbmatlab);
            addpath(rbmatlab);
            curdir = pwd;
            startup_rbmatlab;
            chdir(curdir);
        end
    end
    
end

