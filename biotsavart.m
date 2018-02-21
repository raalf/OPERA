clc
clear

syms xsi eta zeta xm ym zm ex ey ez 
syms A1 A2 B1 B2 C2 C3

ex = [1 0 0];
ey = [0 1 0];
ez = [0 0 1];

r = xm*ex + ym*ey + zm*ez;
s = xsi*ex + eta*ey + zeta*ez;

Gamma = 0.5*A1*eta^2 + A2*eta + 0.5*B1*xsi^2 + B2*xsi + C2*xsi*eta + C3;
gamma = gradient(Gamma, [xsi, eta, zeta]);

zeta = 0;
tr = (cross((r - s), gamma))./(abs(r - s).^3);

syms eta_le eta_te
tr1 = int(tr, eta, eta_le, eta_te);

pretty(tr1)
