clc
clear

syms x_m y_m z_m xi eta zeta A1 A2 B1 B2 C2 C3 J1 J2 J3 J4 J5 J6 xi_1 ...
    xi_2 xi_3 eta_1 eta_2 eta_3 real

r = [x_m, y_m, z_m];
s = [xi, eta, 0];

Gamma = 0.5*A1*(eta^2) + A2*eta + 0.5*B1*(xi^2) + B2*xi + C2*xi*eta + C3;
omega = gradient(Gamma, [xi, eta, zeta]);
omega = [omega(2), -omega(1), omega(3)];

C = (eta_3 - eta_2)/(xi_3 - xi_2);
D_LE = eta_2 - ((xi_2*(eta_3 - eta_2))/(xi_3 - xi_2));
eta_LE = C*xi + D_LE;

E = (eta_3 - eta_1)/(xi_3 - xi_1);
D_TE = eta_1 - ((xi_1*(eta_3 - eta_1))/(xi_3 - xi_1));
eta_TE = E*xi + D_TE;

V = cross(r-s, omega)/(norm(r-s,2)^3);
q = (norm(r-s,2)^3);

u = z_m*(B1*J2 + C2*J4 + B2*J1);
v = z_m*(A1*J4 + C2*J2 + A2*J1);
w = A1*J5 - A1*y_m*J4 - C2*x_m*J4 + 2*C2*J6 + A2*J4 + B1*J3 - B1*x_m*J2 - C2*y_m*J2 + B2*J2 - A2*y_m*J1 -  B2*x_m*J1;

D_kin = equationsToMatrix([u;v;w], [A1; A2; B1; B2; C2; C3]);

assume(xi_1 < xi_3);
assume(xi_1 == xi_2);
assume(eta_2 > eta_1);

