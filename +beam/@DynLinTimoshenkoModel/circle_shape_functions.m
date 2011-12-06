function N = circle_shape_functions(this, R, s, B)

% Wertet die Basisfunktionen für ein Kreiselement aus:
% Auf dem Element werden konstante Strains angenommen => DGLs der Form
%   x' = Ax + c2
% müssen gelöst werden. x = (u1, u2, u3, theta1, theta2, theta3). Es ergibt
% sich eine Lösung der Form:
%   x(s) = C(s) * c. Wobei C 6x12 und c ein Vektor mit Konstanten der Länge 12
% ist. (c = [x(0); c2])
% Um nun Basisfunktionen zu generieren müssen die Lösungen x(s) mit den
% Werten an den beiden Enden des Elements verbunden werden:
% v = [x(0); x(L)] = [C(0); C(L)] * c =: inv(B) * c

%       => x(s) = C(s)*B * v =: N(s) * [x(0); x(L)] !

% Da M für jedes Element konstant ist (hängt nur von L ab!), könnte M im
% Voraus berechnet werden.

% Anmerkung: Die Basisfunktionen sind so gebaut, dass die Ableitung
% (d/ds u = u' + Gu, etc) eben jene Konstanten sind, die in der zweiten
% Hälfte von c stehen (c2 s.o.): d/ds N(s)*v = (B*v)(7:12)
%       (d/ds x(s) = d/ds (C(s))*c = c2 ?)


% C = @(s) [
%     cos(s/R)    sin(s/R)    0   0               0           R-R*cos(s/R)    R*sin(s/R)      R-R*cos(s/R)    0   0                   0                   R*s-R^2*sin(s/R);
%     -sin(s/R)   cos(s/R)    0   0               0           R*sin(s/R)      R*cos(s/R)-R    R*sin(s/R)      0   0                   0                   R^2-R^2*cos(s/R);
%     0           0           1   R-R*cos(s/R)    -R*sin(s/R) 0               0               0               s   R*s-R^2*sin(s/R)    R^2*cos(s/R)-R^2    0;
%     0           0           0   cos(s/R)        sin(s/R)    0               0               0               0   R*sin(s/R)          R-R*cos(s/R)        0;
%     0           0           0   -sin(s/R)       cos(s/R)    0               0               0               0   R*cos(s/R)-R        R*sin(s/R)          0;
%     0           0           0   0               0           1               0               0               0   0                   0                   s];

% B = inv([C(R, 0); C(R, s_end)]);

co = cos(s/R);
si = sin(s/R);
RmRco = R-R*co;
R2mR2co = R*RmRco;
Rsi = R*si;
RsmR2si = R*s - R^2*si;
Cs =  [
    co	si	0	0       0       RmRco	Rsi     RmRco	0	0           0           RsmR2si;
    -si	co	0   0       0       Rsi     -RmRco	Rsi     0   0           0           R2mR2co;
    0	0	1   RmRco	-Rsi	0       0       0       s   RsmR2si     -R2mR2co    0;
    0	0	0   co      si      0       0       0       0   Rsi         RmRco       0;
    0	0	0   -si     co      0       0       0       0   -RmRco      Rsi         0;
    0	0	0   0       0       1       0       0       0   0           0           s];
% Cs =  [
%     co	si	0	0       0       R-R*co	R*si	R-R*co	0	0           0           R*s-R^2*si;
%     -si	co	0   0       0       R*si	R*co-R	R*si	0   0           0           R^2-R^2*co;
%     0	0	1   R-R*co	-R*si	0       0       0       s   R*s-R^2*si	R^2*co-R^2	0;
%     0	0	0   co      si      0       0       0       0   R*si        R-R*co      0;
%     0	0	0   -si     co      0       0       0       0   R*co-R      R*si        0;
%     0	0	0   0       0       1       0       0       0   0           0           s];
N = Cs*B;
