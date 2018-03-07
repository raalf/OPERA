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
tr = (cross((r - s), gamma))./(abs(r - s).^3);

%% Leading and trailing edge integration limits
% Eta location of leading and trailing edge changes linearly with xi
syms xi_1 xi_2 xi_3 eta_1 eta_2 eta_3 real
eta_le = eta_2 + (xi - xi_2)*((eta_3 - eta_2)/(xi_3 - xi_2));
eta_te = eta_1 + (xi - xi_1)*((eta_3 - eta_1)/(xi_3 - xi_1));

assumeAlso(xi_3 > xi_2);
assumeAlso(xi_2 == xi_1);
assumeAlso(eta_2 > eta_1);

%% Integrating Biot-Savart from trailing to leading edge
tr1(1) = int(tr(1), eta, eta_te, eta_le);
tr1(2) = int(tr(2), eta, eta_te, eta_le);
tr1(3) = int(tr(3), eta, eta_te, eta_le);

%% Integrating from left to right (xi_1 to xi_3)

% % Xi direction w_ind
% kids = children(tr1(1));
% for i = 1:length(kids)
% %     clc
%     disp([num2str(i) '/' num2str(length(kids))])
%     tr21_pre(i) = int(simplify(kids(i)), xi, xi_1, xi_3);
% %     pretty(tr21_pre(i))
% end
% % Summing the child integrals back up and evaluating from TE to LE
% tr2(1) = sum(tr21_pre);

% % Eta direction w_ind
% kids = children(expand(tr1(2)));
% for i = 1:length(kids)
%     disp(i)
%     tr22_pre(i) = int(simplify(kids(i)), xi, xi_1, xi_3);
% end
% % The sign of the children terms are based on the sign of xi - xi_p
% tr22_pre = tr22_pre.*sign(xi - xi_p); 
% % Summing the child integrals back up and evaluating from TE to LE
% tr22_comb = sum(tr22_pre);
% tr2(2) = subs(tr22_comb, xi, xi_3) - subs(tr22_comb, xi, xi_1);

% Zeta direction w_ind (easy)
tr2(3) = int(tr1(3), xi, xi_1, xi_3);

matlabFunction(tr2(3),'File','w_zeta')

% fid = fopen('bs2.txt', 'wt');
% for i = 1:2
%     
%     kids = children(tr2(i));
%     
%     for jj = 1:length(kids)
%         if jj < length(kids)
%             fprintf(fid, '+ (%s) ... \n', char(kids(jj)));
%         elseif jj == length(kids)
%             fprintf(fid, '+ (%s); \n', char(kids(jj))); 
%         end
%     end
%     
%     fprintf(fid, '\n\n\n\n');
%     
% end
% 
% fprintf(fid, '%s; \n\n', char(tr2(3)));
% fclose(fid);













