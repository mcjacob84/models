classdef StructureElement < handle
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
        PointsIdx;
        
        % Material vector
        Material;
    end
    
    properties(SetAccess = protected)
        % Indizes (lokal pro Knoten) in die die lokalen Matrizen assembliert werden
        % Konvention: Freiheitsgrade pro Knoten (u1, u2, u3, phi1, phi2, phi3, T, [...])
%         MatrixDofIndices;
                
        % Länge des Elements
        Length;
        
        % The model that contains the structure element
        Model;
        
        % Transformationsmatrix für das lokale Koordinatensystem
        T;
        
        % Effektive Stoffkonstanten
        c;
        
        c_theta = [];
        kappa = [];
        alphaA = [];
    end
    
    properties(Access=protected)
        % The global transformation matrix (?)
        %
        % Is assigned in the subclasses' initialize method.
        TG;
    end
    
    methods
        function this = StructureElement(model, material, pointsidx)
            this.Model = model;
            this.PointsIdx = pointsidx;
            this.Material = material;
        end
    end
    
    methods(Abstract)
        M = getLocalMassMatrix(this);
        K = getLocalStiffnessMatrix(this);
        f = getLocalForce(this, gravity);
        [K, R, U_pot] = getLocalTangentials(this, u);
        globIdx = getGlobalIndices(this);
    end
    
end
