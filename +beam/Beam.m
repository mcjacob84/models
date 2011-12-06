classdef Beam < handle
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
    
    properties(Constant)
        % Außendurchmesser (m)
        ROHR_d_a = 457e-3;
        % Wandstärke (m)
        ROHR_s = 40e-3;
        % Isolierungsdicke (m)
        ROHR_iso = 400e-3;
        % Manteldicke (m)
        ROHR_mantel = 1e-3;
        % Dichte des Stahls (kg/m³)(eingelesen!)
        ROHR_rho = c(1, 1);
        % Dichte der Isolierung (kg/m³)
        ROHR_rho_iso = 100;
        % Dichte des Mantels (kg/m³)
        ROHR_rho_mantel = 7850;
        % Dichte des Mediums (kg/m³)
        ROHR_rho_med = 20;
        % Querkontraktionszahl
        ROHR_ny = 0.3;
        % E-Modul (N/m²)(eingelesen!)
        ROHR_E = c(1, 3);
        % Rohrinnendruck (N/m²)
        ROHR_p = 50 * 1e5;
    end
    
    methods
    end
    
end