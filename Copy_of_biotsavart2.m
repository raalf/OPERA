clc
clear

syms xsi eta zeta xm ym zm ex ey ez 
syms A1 A2 B1 B2 C2 C3

ex = [1 0 0];
ey = [0 1 0];
ez = [0 0 1];

zeta = 0;
r = xm*ex + ym*ey + zm*ez;
s = xsi*ex + eta*ey + zeta*ez;

% Gamma = 0.5*B1*xsi^2 + B2*xsi + C3;
Gamma = 0.5*A1*eta^2 + A2*eta + 0.5*B1*xsi^2 + B2*xsi + C2*xsi*eta + C3;
gamma = [gradient(Gamma, [xsi, eta]); 0];
gamma = [-gamma(2); gamma(1); gamma(3)]; % Shed vorticity is perpendicular?

tr = (cross((r - s), gamma))./(abs(r - s).^3);

syms xsi_1 xsi_2 xsi_3 eta_1 eta_2 eta_3 xsi_12

eta_le = eta_2 + (xsi - xsi_2)*((eta_3 - eta_2)/(xsi_3 - xsi_2));
eta_te = eta_1 + (xsi - xsi_1)*((eta_3 - eta_1)/(xsi_3 - xsi_1));

tr1(1) = int(tr(1), eta, eta_te, eta_le);
tr1(3) = int(tr(3), eta, eta_te, eta_le);
% CHECK THIS vvvvvvvvvvvvvv
tr2_pre = -(A2*eta*zm)/(abs(eta - ym)^3) - (eta*A1*eta*zm)/(abs(eta - ym)^3) - (C2*eta*xsi*zm)/(abs(eta - ym)^3);
tr1(2) = subs(tr2_pre, eta, eta_le) - subs(tr2_pre, eta, eta_te);

tr11_c = children(tr1(1));
tr12_c = children(tr1(2));
tr13_c = children(tr1(3));

for i = 1:3
    kids = children(tr1(i));
    for j = 1:size(kids,2)
       disp(['Direction ', num2str(i), ' Term ', num2str(j)]) 
       t(i,j) = int(kids(j), xsi);
    end
    kids = [];
end
