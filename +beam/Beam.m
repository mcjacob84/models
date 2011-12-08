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
    
    properties(Constant)

    end
    
    methods
         function this = Beam(material, points)
            this = this@models.beam.StructureElement(material, points);
            
                % Materialsatz
    RO(i).Mat = this.RO_raw(i,1);
    
    % Anfangs- und Endpunkt
    RO(i).p = [this.RO_raw(i,2) this.RO_raw(i,3)];
    
    % Längenberechnung
    dx = (this.Points(RO(i).p(2), 1) - this.Points(RO(i).p(1), 1));
    dy = (this.Points(RO(i).p(2), 2) - this.Points(RO(i).p(1), 2)); 
    dz = (this.Points(RO(i).p(2), 3) - this.Points(RO(i).p(1), 3)); 
    RO(i).l = sqrt( dx^2 + dy^2 + dz^2 );
    
    % Lokales Koordinatensystem
    e_x = [dx / RO(i).l; dy / RO(i).l; dz / RO(i).l];
    e_y = [-e_x(2) e_x(1) 0]';
    if norm(e_y) == 0
        e_y = [0 e_x(3) -e_x(2)]';
    end
    e_z = [e_x(2)*e_y(3) - e_x(3)*e_y(2);
           e_x(3)*e_y(1) - e_x(1)*e_y(3);
           e_x(1)*e_y(2) - e_x(2)*e_y(1)];
    e_y = e_y / norm(e_y);
    e_z = e_z / norm(e_z);
    RO(i).T = [e_x e_y e_z];
    
    % Effektive Konstanten
%   <rho>	<A>     <E>     <Iy/In>     <Iz/Ib>     <It>        <G>     <k>	<c_th>	<kappa>	<alpha>
%   1       2       3       4           5           6           7       8   9       10      11
    mat = RO(i).Mat;
    RO(i).c(1) = this.mat(mat,3) * this.mat(mat,4);                           % c1            = E*I
    RO(i).c(2) = RO(i).c(1)^2;                                  % c2 = c1^2     = (E*I)^2
    RO(i).c(3) = this.mat(mat,8)*this.mat(mat,7)*this.mat(mat,2)*RO(i).l;            % c3            = G*As*L
    RO(i).c(4) = RO(i).c(3)^2;                                  % c4 = c3^2     = (G*As*L)^2
    RO(i).c(5) = RO(i).c(3) * RO(i).l;                          % c5 = c3*L     = G*As*L^2
    RO(i).c(6) = RO(i).c(5)^2;                                  % c6 = c5^2     = (G*As*L^2)^2
    RO(i).c(7) = RO(i).c(1) * RO(i).c(5);                       % c7 = c1*c5    = E*I*G*As*L^2
    RO(i).c(8) = RO(i).c(1) * this.mat(mat,8)*this.mat(mat,7)*this.mat(mat,2);       % c8 = c1*G*As  = E*I*G*As
    RO(i).c(9) = this.mat(mat,1) * this.mat(mat,2);                           % c9 = rho*A
    RO(i).c(10) = this.mat(mat,1) * this.mat(mat,2) * 9.81;                   %  q = rho*A*Ortsfaktor
    RO(i).c(11) = this.mat(mat,1) * this.mat(mat,4);                          % c11= rho*I
    RO(i).c(12) = RO(i).c(9) * RO(i).l / 6;                     % c12= c9*L/6   = rho*A*L/6
    RO(i).c(13) = this.mat(mat,3) * this.mat(mat,2) / RO(i).l;                % c13= E*A/L
    RO(i).c(14) = this.mat(mat,7)*this.mat(mat,6) / RO(i).l;                  % c14= G*It/L
    RO(i).c(15) = this.mat(mat,1)*this.mat(mat,6)*RO(i).l / 6;                % c15= rho*It*L/6
    
    RO(i).q_lok = (RO(i).c(9) + ROHR_q_plus) * (RO(i).T' * gravity);
    RO(i).c_theta = this.mat(mat,9);
    RO(i).kappa = this.mat(mat,10);
    RO(i).alphaA = this.mat(mat,11) * this.mat(mat,2) * this.mat(mat,3);             % alpha*A*E

        end
    end
    
end