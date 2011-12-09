classdef Material < handle
% Material: 
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
    
    properties
        % Außendurchmesser (m)
        d_a = 457e-3;
        % Wandstärke (m)
        s = 40e-3;
        % Isolierungsdicke (m)
        iso = 400e-3;
        % Manteldicke (m)
        mantel = 1e-3;
        % Dichte des Stahls (kg/m³) (Konstruktor!)
        rho;
        % Dichte der Isolierung (kg/m³)
        rho_iso = 100;
        % Dichte des Mantels (kg/m³)
        rho_mantel = 7850;
        % Dichte des Mediums (kg/m³)
        rho_med = 20;
        % Querkontraktionszahl
        ny = 0.3;
        % E-Modul (N/m²) (Konstruktor!)
        E;
        % Rohrinnendruck (N/m²)
        p = 50 * 1e5;
    end
    
    properties(SetAccess = protected)
        A;
        Iy;
        Iz;
        It;
        G;
        k;
        q_plus;
    end
    
    methods
        function this = Material(rho, E)
            % @todo q_pkus in effektive dichte umrechnen
            this.rho = rho;
            this.E = E;
            
            % Innen-/Außendurchmesser des Rohrs
            r_a = 0.5 * this.d_a;
            r_i = r_a - this.s;
            % Querschnittsfläche für Balken
            this.A = pi * ( r_a^2 - r_i^2 );
            % Flächenträgheitsmoment für Balken
            this.Iy = 0.25 * pi * ( r_a^4 - r_i^4 );
            this.Iz = this.Iy;
            
            % Torsionsträgheitsmoment für Balken
            this.It = 2 * Iy;
            % Schubmodul für Balken
            this.G = E / ( 2*(1+this.ny) );
            % Schubkorrekturfaktor für Balken
            m_tmp = r_i / r_a;
            this.k = 6*(1+this.ny)*(1+m_tmp^2)^2 / ( (7+6*this.ny)*(1+m_tmp^2)^2 + (20+12*this.ny)*m_tmp^2);

            % Berechnung der durch Medium und Dämmung verursachten zusätzlichen Steckenlast
            this.q_plus = pi * ( r_i^2 * this.rho_med + ( (r_a + this.iso)^2 - r_a^2 ) * this.rho_iso + ( (r_a + this.iso + this.mantel)^2 - (r_a + this.iso)^2) * this.rho_mantel );
        end
        
    end
    
end