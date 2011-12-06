function N = beam_shape_functions_derivative(this, s, L, c)

% Wertet (analytische) Basisfunktion für einen geraden Balken der Länge L und den Stoffkonstanten c an der Stelle s aus.
% c1            = E*I
% c2 = c1^2     = (E*I)^2
% c3            = G*As*L
% c4 = c3^2     = (G*As*L)^2
% c5 = c3*L     = G*As*L^2
% c6 = c5^2     = (G*As*L^2)^2
% c7 = c1*c5    = E*I*G*As*L^2
% c8 = c1*G*As  = E*I*G*As
% c9 = rho*A
% c10= rho*A*Ortsfaktor
% c11= rho*I
% c12= c9*L/6   = rho*A*L/6
% c13= E*A/L
% c14= G*It/L
% c15= rho*It*L/6

GAs = c(3)/L;
GAsL = c(3);
GAsLL = c(5);
EIy = c(1);
EIz = c(1);

N_u =  [-1/L 0 0 0 0 0 1/L 0 0 0 0 0];
N_phi = [0 0 0 -1/L 0 0 0 0 0 1/L 0 0];

nenner_y = 1/(L * (12*EIy + GAsLL));
nenner_z = 1/(L * (12*EIz + GAsLL));
nenner2_y = 1/(12*EIy + GAsLL);
nenner2_z = 1/(12*EIz + GAsLL);

V_v1 = nenner_y * (6*GAs * s^2 - 6*GAsL * s - 12*EIy);
V_theta1 = nenner_y * (3*GAsL * s^2 - 4*s*(3*EIy + GAsLL) + L*(6*EIy + GAsLL));
V_v2 = -nenner_y * (6*GAs * s^2 - 6*GAsL * s - 12*EIy);
V_theta2 = nenner_y * (3*GAsL * s^2 + 2*s*(6*EIy - GAsLL) - 6*EIy*L);

N_v = [0 V_v1 0 0 0 V_theta1 0 V_v2 0 0 0 V_theta2];


W_w1 = nenner_z * (6*GAs * s^2 - 6*GAsL * s - 12*EIz);
W_psi1 = -nenner_z * (3*GAsL * s^2 - 4*s*(3*EIz+GAsLL) + L*(6*EIz + GAsLL));
W_w2 = -nenner_z * (6*GAs * s^2 - 6*GAsL * s - 12*EIz);
W_psi2 = -nenner_z * (3*GAsL * s^2 + 2*s*(6*EIz - GAsLL) - 6*EIz*L);

N_w = [0 0 W_w1 0 W_psi1 0 0 0 W_w2 0 W_psi2 0];

Psi_w1   = nenner_z *6*GAs * (L - 2*s);
Psi_psi1 = nenner2_z * 6*GAs * s - nenner_z * 4*(3*EIz + GAsLL);
Psi_w2   = nenner_z * 6*GAs * (2*s - L);
Psi_psi2 = nenner2_z * 6*GAs * s + nenner_z * 2*(6*EIz - GAsLL);

N_psi = [0 0 Psi_w1 0 Psi_psi1 0 0 0 Psi_w2 0 Psi_psi2 0];

Theta_v1 = nenner_y * 6*GAs * (2*s - L);
Theta_theta1 = nenner2_y * 6*GAs * s - nenner_y * 4*(3*EIy + GAsLL);
Theta_v2 = nenner_y * 6*GAs * (L - 2*s);
Theta_theta2 = nenner2_y * 6*GAs * s + nenner_y * 2*(6*EIy - GAsLL);

N_theta = [0 Theta_v1 0 0 0 Theta_theta1 0 Theta_v2 0 0 0 Theta_theta2];

N = [N_u; N_v; N_w; N_phi; N_psi; N_theta];
