classdef Beam < models.beam.StructureElement
% Beam: 
%
%
%
% @author Daniel Wirtz @date 2011-12-05
%
% @new{0,6,dw,2011-12-05} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    properties
        % Speichert den anteil des jeweiligen Elements an der Gesamtsystemlänge
        %
        % Siehe preprocess_data
        split;
    end
    
    methods
        function this = Beam(model, material, pointsidx)
            this = this@models.beam.StructureElement(model, material, pointsidx);
        end
        
        function globIdx = getGlobalIndices(this)
            % Returns the global indices in the vector of DoFs for this
            % element
            %
            % Identical for straight and curved beams.
            data = this.Model.data;
            p = this.PointsIdx;
            index_0_glob = 7*data.knot_index(p(1))-6 : 7*data.knot_index(p(1))-1;
            index_l_glob = 7*data.knot_index(p(2))-6 : 7*data.knot_index(p(2))-1;
            globIdx = [index_0_glob index_l_glob];
        end
    end
    
    methods(Access=protected)
        function initialize(this)
%             this.c_theta = this.Material(9);
%             this.kappa = this.Material(10);
%             this.alphaA = this.Material(11) * this.Material(2) * this.Material(3); % alpha*A*E
        end
    end
    
end