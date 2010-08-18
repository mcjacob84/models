classdef RBMDataContainer < handle
    %RBMDATACONTAINER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        RBMData;
    end
    
    methods
        function this = RBMDataContainer(rbmatlabdata)
            this.RBMData = rbmatlabdata;
        end
    end
    
end

