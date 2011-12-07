classdef StructureElement < 
% StructureElement: 
%
%
%
% @author Christoph Strohmeyer @date 2011-12-06
%
% @new{0,6,CS,2011-12-06} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetAccess = private)
        % Punkt-Index-Array
        Points;
        
        % Material
        Material;
    end
    
    properties(SetAccess = protected)
    
        % Indizes (lokal pro Knoten) in die die lokalen Matrizen assembliert werden
        % Konvention: Freiheitsgrade pro Knoten (u1, u2, u3, phi1, phi2, phi3, T, [...])
        MatrixDofIndices;
                
        % Länge des Elements
        Length;
    end
    
    properties(Access = protected)
        % Effektive Stoffkonstanten
        c;
        
        % Transformationsmatrix für das lokale Koordinatensystem
        T;
    end
    
    methods
        function this = StructureElement(material, points)
            this.Points = points;
            this.Material = material;
        end
    end
    
    methods(Abstract)
        [M, K, f] = getLocalMatrices(this);
    end
    
end
