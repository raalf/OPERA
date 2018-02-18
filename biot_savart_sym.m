clc
clear

syms u Gamma gamma x x_p xsi eta zeta xsi_p eta_p zeta_p

x_p = [xsi_p eta_p 0];
x = [xsi eta zeta];

syms A1 A2 B1 B2 C2 C3

Gamma = 0.5*A1*xsi_p^2 + A2*xsi_p + 0.5*B2*eta_p^2 + B2*eta_p + C2*eta_p*xsi_p + C3;
gamma = gradient(Gamma, [xsi_p, eta_p, zeta_p]);

r = x - x_p;
u = int(int( cross(gamma, r) ./ abs(r).^3, xsi_p), eta_p)