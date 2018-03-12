clc
clear

reset(symengine)

syms xi eta zeta xi_p eta_p zeta_p ex ey ez real
syms A1 A2 B1 B2 C2 C3 real

%% Local reference frame axis (xi, eta, zeta)
ex = [1 0 0];
ey = [0 1 0];
ez = [0 0 1];

zeta = 0; % For circulation and vorticity, zeta is always 0

% Vector from field point to surface element
r = xi_p*ex + eta_p*ey + zeta_p*ez;
% Vector from element origin to surface element
s = xi*ex + eta*ey + zeta*ez;

%% Circulation and vorticity
Gamma = 0.5*A1*eta^2 + A2*eta + 0.5*B1*xi^2 + B2*xi + C2*xi*eta + C3;
gamma = [gradient(Gamma, [xi, eta]); 0];
gamma = [-gamma(2); gamma(1); gamma(3)]; % Shed vorticity is perpendicular

%% Biot-Savart integral
% According to Green's theorem
tr = (cross((r - s), gamma))./(norm(r - s).^3);

%% Leading and trailing edge integration limits
% Eta location of leading and trailing edge changes linearly with xi
syms xi_1 xi_2 xi_3 eta_1 eta_2 eta_3 real
eta_le = eta_2 + (xi - xi_2)*((eta_3 - eta_2)/(xi_3 - xi_2));
eta_te = eta_1 + (xi - xi_1)*((eta_3 - eta_1)/(xi_3 - xi_1));

assumeAlso(xi_3 > xi_2);
assumeAlso(xi_2 == xi_1);
assumeAlso(eta_2 > eta_1);

coeff = [A1; A2; B1; B2; C2; C3];

D_kin = equationsToMatrix(tr, coeff)

%% Integrating Biot-Savart from trailing to leading edge

% fid = fopen('bs1.txt', 'wt');
% fprintf(fid, 'Xi term \n\n');
% 
% % Xi term
% kids = children(expand(tr(1)));
% for i = 1:length(kids)
%     clc
%     disp(['Xi-dir ', num2str(i) '/' num2str(length(kids))])
%     tr11_pre(i) = int(simplify(kids(i)), eta, eta_te, eta_le);
%     fprintf(fid, '%s \n', char(tr11_pre(i)));
%     pretty(tr11_pre(i))
% end
% fprintf(fid, '\n');
% % Summing the child integrals back up and evaluating from TE to LE
% % tr1(1) = sum(tr11_pre);
% % fprintf(fid, '%s \n\n', char(tr1(1)));
% 
% fprintf(fid, 'Eta term \n\n');
% % Eta term
% kids = children(expand(tr(2)));
% for i = 1:length(kids)
%     clc
%     disp(['Eta-dir ', num2str(i) '/' num2str(length(kids))])
%     tr12_pre(i) = int(simplify(kids(i)), eta, eta_te, eta_le);
%     fprintf(fid, '%s \n', char(tr12_pre(i)));
%     pretty(tr12_pre(i))
% end
% fprintf(fid, '\n');
% % Summing the child integrals back up and evaluating from TE to LE
% % tr1(2) = sum(tr12_pre);
% % fprintf(fid, '%s \n\n', char(tr1(2)));
% 
% % Zeta term
% kids = children(expand(tr(3)));
% for i = 1:length(kids)
%     clc
%     disp(['Zeta-dir ', num2str(i) '/' num2str(length(kids))])
%     tr13_pre(i) = int(simplify(kids(i)), eta, eta_te, eta_le);
%     fprintf(fid, '%s \n', char(tr13_pre(i)));
%     pretty(tr13_pre(i))
% end
% fprintf(fid, '\n');
% % Summing the child integrals back up and evaluating from TE to LE
% % tr1(1) = sum(tr11_pre);
% % fprintf(fid, '%s \n\n', char(tr1(1)));
% 
% fclose(fid);
% 
% %% Integrating from left to right (xi_1 to xi_3)












