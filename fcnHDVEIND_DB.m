function [infl_loc] = fcnHDVEIND_DB(dvenum, dvetype, fpg, matPLEX, matROTANG, matCONTROL, vecBI, ztol)
warning('on')
tol = 1e-2;

del_h = 0.01;

fpl = fcnGLOBSTAR(fpg - matCONTROL(dvenum,:), matROTANG(dvenum,:));
len = size(fpl,1);

x_m = fpl(:,1);
y_m = fpl(:,2);
h = fpl(:,3);

xi_1 = permute(matPLEX(1,1,dvenum),[3 2 1]);
xi_2 = permute(matPLEX(2,1,dvenum),[3 2 1]);
xi_3 = permute(matPLEX(3,1,dvenum),[3 2 1]);

eta_1 = permute(matPLEX(1,2,dvenum),[3 2 1]);
eta_2 = permute(matPLEX(2,2,dvenum),[3 2 1]);
eta_3 = permute(matPLEX(3,2,dvenum),[3 2 1]);

%%
% Checking which elements are on the element
C = (eta_3 - eta_2)./(xi_3 - xi_2);
D_LE = eta_2 - ((xi_2.*(eta_3 - eta_2))./(xi_3 - xi_2));
E = (eta_3 - eta_1)./(xi_3 - xi_1);
D_TE = eta_1 - ((xi_1.*(eta_3 - eta_1))./(xi_3 - xi_1));

le_eta = C.*x_m + D_LE;
te_eta = E.*x_m + D_TE;
xi_left = min([xi_1, xi_3],[],2);
xi_right = max([xi_1, xi_3],[],2);

margin = 1e-3;
idx_on_element = y_m >= te_eta - margin & y_m <= le_eta + margin & x_m >= xi_left - margin & x_m <= xi_right + margin & abs(h) <= margin;

%%
% hFig1 = figure(1);
% clf(1);
% patch([xi_1(1) xi_2(1) xi_3(1)], [eta_1(1) eta_2(1) eta_3(1)],'w');
% grid minor
% box on
% axis tight
% axis equal
% xlabel('Xi-Dir');
% ylabel('Eta-Dir');
% zlabel('Zi-Dir');
% for i = 1:size(fpl,1)
%     hold on
%     scatter3(fpl(i,1), fpl(i,2), fpl(i,3), 'rs')
%     hold off
% end

% Call, xi-eta, edge (1, 2, or 3), point (1 or 2)
epts = nan(len,2,3,2);
epts(:,:,1,1:2) = cat(4, [xi_1 eta_1], [xi_2 eta_2]);
epts(:,:,2,1:2) = cat(4, [xi_2 eta_2], [xi_3 eta_3]);
epts(:,:,3,1:2) = cat(4, [xi_3 eta_3], [xi_1 eta_1]);

% Edge direction and normal (d and n) (depth is edge number)
nu_d = epts(:,:,:,2) - epts(:,:,:,1);
nu_d = nu_d./sqrt(sum(nu_d.^2,2));

nu_n = nu_d(:,[2 1],:);
nu_n(:,1,:) = nu_n(:,1,:).*-1;

% Calculating parameters
L1 = dot(nu_d, fpl(:,1:2) - epts(:,:,1:3,1), 2);
L2 = dot(nu_d, fpl(:,1:2) - epts(:,:,1:3,2), 2);
a = dot(nu_n, mean(epts,4) - fpl(:,1:2), 2);

d_H = min(mean([sqrt(L1.^2 + a.^2) sqrt(L2.^2 + a.^2)],2),[],3);

g = sqrt(a.^2 + h.^2);
c1 = g.^2 + abs(h).*sqrt(L1.^2 + g.^2);
c2 = g.^2 + abs(h).*sqrt(L2.^2 + g.^2);

%% E
% MXK = 5, MXQ = 5, MXFK = 16 + 5 - 2 = 19;

% PAGE 183-184
rho1 = sqrt(L1.^2 + g.^2);
rho2 = sqrt(L2.^2 + g.^2);

% a.)
x2 = rho2.^-2;
x1 = rho1.^-2;
A2 = (epts(:,1,:,2) - x_m)./rho2;
A1 = (epts(:,1,:,1) - x_m)./rho1;

E211 = A2 - A1; 
E213 = A2.*x2 - A1.*x1; 
E215 = (x2 + x1).*E213 - (x1.*x2.*E211); 
E217 = (x2 + x1).*E215 - (x1.*x2.*E213);
E219 = (x2 + x1).*E217 - (x1.*x2.*E215);
E2111 = (x2 + x1).*E219 - (x1.*x2.*E217);
E2113 = (x2 + x1).*E2111 - (x1.*x2.*E219);
E2115 = (x2 + x1).*E2113 - (x1.*x2.*E2111);
E2117 = (x2 + x1).*E2115 - (x1.*x2.*E2113);
E2119 = (x2 + x1).*E2117 - (x1.*x2.*E2115);
E2121 = (x2 + x1).*E2119 - (x1.*x2.*E2117);
E2123 = (x2 + x1).*E2121 - (x1.*x2.*E2119);
E2125 = (x2 + x1).*E2123 - (x1.*x2.*E2121);
E2127 = (x2 + x1).*E2125 - (x1.*x2.*E2123);
E2129 = (x2 + x1).*E2127 - (x1.*x2.*E2125);
E2131 = (x2 + x1).*E2129 - (x1.*x2.*E2127);
E2133 = (x2 + x1).*E2131 - (x1.*x2.*E2129);

% I = 3;
% E215 = A2.*(x2.^(I - 1)) - A1.*(x1.^(I - 1)); 
% I = 4;
% E217 = A2.*(x2.^(I - 1)) - A1.*(x1.^(I - 1)); 
% I = 5;
% E219 = A2.*(x2.^(I - 1)) - A1.*(x1.^(I - 1)); 
% I = 6;
% E2111 = A2.*(x2.^(I - 1)) - A1.*(x1.^(I - 1)); 
% I = 7;
% E2113 = A2.*(x2.^(I - 1)) - A1.*(x1.^(I - 1)); 
% I = 8;
% E2115 = A2.*(x2.^(I - 1)) - A1.*(x1.^(I - 1)); 
% I = 9;
% E2117 = (x2 + x1).*E2115 - (x1.*x2.*E2113);

% b.)
A2 = (epts(:,2,:,2) - y_m)./rho2;
A1 = (epts(:,2,:,1) - y_m)./rho1;

E121 = A2 - A1; 
E123 = A2.*x2 - A1.*x1; 
E125 = (x2 + x1).*E123 - (x1.*x2.*E121); 
E127 = (x2 + x1).*E125 - (x1.*x2.*E123);
E129 = (x2 + x1).*E127 - (x1.*x2.*E125);
E1211 = (x2 + x1).*E129 - (x1.*x2.*E127);
E1213 = (x2 + x1).*E1211 - (x1.*x2.*E129);
E1215 = (x2 + x1).*E1213 - (x1.*x2.*E1211);
E1217 = (x2 + x1).*E1215 - (x1.*x2.*E1213);
E1219 = (x2 + x1).*E1217 - (x1.*x2.*E1215);
E1221 = (x2 + x1).*E1219 - (x1.*x2.*E1217);
E1223 = (x2 + x1).*E1221 - (x1.*x2.*E1219);
E1225 = (x2 + x1).*E1223 - (x1.*x2.*E1221);
E1227 = (x2 + x1).*E1225 - (x1.*x2.*E1223);
E1229 = (x2 + x1).*E1227 - (x1.*x2.*E1225);
E1231 = (x2 + x1).*E1229 - (x1.*x2.*E1227);
E1233 = (x2 + x1).*E1231 - (x1.*x2.*E1229);

% c.)
E11n1 = rho2 - rho1; % I = 1
E12n1 = rho2.*(epts(:,2,:,2) - y_m) - rho1.*(epts(:,2,:,1) - y_m); % I = 2
E13n1 = rho2.*(epts(:,2,:,2) - y_m).^2 - rho1.*(epts(:,2,:,1) - y_m).^2; % I = 3
E14n1 = rho2.*(epts(:,2,:,2) - y_m).^3 - rho1.*(epts(:,2,:,1) - y_m).^3; % I = 4

% d.)
E21n1 = rho2.*(epts(:,1,:,2) - x_m) - rho1.*(epts(:,1,:,1) - x_m); % I = 2
E31n1 = rho2.*(epts(:,1,:,2) - x_m).^2 - rho1.*(epts(:,1,:,1) - x_m).^2; % I = 3
E41n1 = rho2.*(epts(:,1,:,2) - x_m).^3 - rho1.*(epts(:,1,:,1) - x_m).^3; % I = 4

%e.)
E111 = (1./rho2) - (1./rho1);

%% F
% MXK = 5, MXQ = 5, MXFK = 16 + 5 - 2 = 19;
% PAGE 130
% 1.)
F111 = zeros(size(a));

idx = L1 >= 0 & L2 >= 0;
F111(idx) = log((rho2(idx) + L2(idx))./(rho1(idx) + L1(idx)));

idx = L1 < 0 & L2 < 0;
F111(idx) = log((rho1(idx) - L1(idx))./(rho2(idx) - L2(idx)));

% idx = L1 < 0 & L2 >= 0;
idx = (L1 < 0 & L2 >= 0) | (L1 >=0 & L2 < 0);
F111(idx) = log(((rho1(idx) - L1(idx)).*(rho2(idx) + L2(idx)))./(g(idx).^2));

idx = [g >= del_h.*d_H, g < del_h.*d_H];
% 2.)
nu_xi = nu_n(:,1,:);
nu_eta = nu_n(:,2,:);

F113 = zeros(size(F111));
F115 = F113;
F117 = F113;
F119 = F113;
F1111 = F113;
F1113 = F113;
F1115 = F113;
F1117 = F113;
F1119 = F113;

F113(idx(:,1,:)) = (1./(g(idx(:,1,:)).^2)).*(-nu_eta(idx(:,1,:)).*E211(idx(:,1,:)) + nu_xi(idx(:,1,:)).*E121(idx(:,1,:)));
K = 5;
F115(idx(:,1,:)) = (1./((g(idx(:,1,:)).^2).*(K - 2))).*((K - 3).*F113(idx(:,1,:)) - nu_eta(idx(:,1,:)).*E213(idx(:,1,:)) + nu_xi(idx(:,1,:)).*E123(idx(:,1,:)));
K = 7;
F117(idx(:,1,:)) = (1./((g(idx(:,1,:)).^2).*(K - 2))).*((K - 3).*F115(idx(:,1,:)) - nu_eta(idx(:,1,:)).*E215(idx(:,1,:)) + nu_xi(idx(:,1,:)).*E125(idx(:,1,:)));
K = 9;
F119(idx(:,1,:)) = (1./((g(idx(:,1,:)).^2).*(K - 2))).*((K - 3).*F117(idx(:,1,:)) - nu_eta(idx(:,1,:)).*E217(idx(:,1,:)) + nu_xi(idx(:,1,:)).*E127(idx(:,1,:)));
K = 11;
F1111(idx(:,1,:)) = (1./((g(idx(:,1,:)).^2).*(K - 2))).*((K - 3).*F119(idx(:,1,:)) - nu_eta(idx(:,1,:)).*E219(idx(:,1,:)) + nu_xi(idx(:,1,:)).*E129(idx(:,1,:)));
K = 13;
F1113(idx(:,1,:)) = (1./((g(idx(:,1,:)).^2).*(K - 2))).*((K - 3).*F1111(idx(:,1,:)) - nu_eta(idx(:,1,:)).*E2111(idx(:,1,:)) + nu_xi(idx(:,1,:)).*E1211(idx(:,1,:)));
K = 15;
F1115(idx(:,1,:)) = (1./((g(idx(:,1,:)).^2).*(K - 2))).*((K - 3).*F1113(idx(:,1,:)) - nu_eta(idx(:,1,:)).*E2113(idx(:,1,:)) + nu_xi(idx(:,1,:)).*E1213(idx(:,1,:)));
K = 17;
F1117(idx(:,1,:)) = (1./((g(idx(:,1,:)).^2).*(K - 2))).*((K - 3).*F1115(idx(:,1,:)) - nu_eta(idx(:,1,:)).*E2115(idx(:,1,:)) + nu_xi(idx(:,1,:)).*E1215(idx(:,1,:)));
K = 19;
F1119(idx(:,1,:)) = (1./((g(idx(:,1,:)).^2).*(K - 2))).*((K - 3).*F1117(idx(:,1,:)) - nu_eta(idx(:,1,:)).*E2117(idx(:,1,:)) + nu_xi(idx(:,1,:)).*E1217(idx(:,1,:)));

% PAGE 132
% NFK = 16, MXFK = 19
F1135 = zeros(size(F111));
F1133 = F1135;
F1131 = F1135;
F1129 = F1135;
F1127 = F1135;
F1125 = F1135;
F1123 = F1135;
F1121 = F1135;

K = 35;
F1133(idx(:,2,:)) = (1./(K - 3)).*((g(idx(:,2,:)).^2).*(K - 2).*F1135(idx(:,2,:)) + nu_eta(idx(:,2,:)).*E2133(idx(:,2,:)) - nu_xi(idx(:,2,:)).*E1233(idx(:,2,:)));
K = 33;
F1131(idx(:,2,:)) = (1./(K - 3)).*((g(idx(:,2,:)).^2).*(K - 2).*F1133(idx(:,2,:)) + nu_eta(idx(:,2,:)).*E2131(idx(:,2,:)) - nu_xi(idx(:,2,:)).*E1231(idx(:,2,:)));
K = 31;
F1129(idx(:,2,:)) = (1./(K - 3)).*((g(idx(:,2,:)).^2).*(K - 2).*F1131(idx(:,2,:)) + nu_eta(idx(:,2,:)).*E2129(idx(:,2,:)) - nu_xi(idx(:,2,:)).*E1229(idx(:,2,:)));
K = 29;
F1127(idx(:,2,:)) = (1./(K - 3)).*((g(idx(:,2,:)).^2).*(K - 2).*F1129(idx(:,2,:)) + nu_eta(idx(:,2,:)).*E2127(idx(:,2,:)) - nu_xi(idx(:,2,:)).*E1227(idx(:,2,:)));
K = 27;
F1125(idx(:,2,:)) = (1./(K - 3)).*((g(idx(:,2,:)).^2).*(K - 2).*F1127(idx(:,2,:)) + nu_eta(idx(:,2,:)).*E2125(idx(:,2,:)) - nu_xi(idx(:,2,:)).*E1225(idx(:,2,:)));
K = 25;
F1123(idx(:,2,:)) = (1./(K - 3)).*((g(idx(:,2,:)).^2).*(K - 2).*F1125(idx(:,2,:)) + nu_eta(idx(:,2,:)).*E2123(idx(:,2,:)) - nu_xi(idx(:,2,:)).*E1223(idx(:,2,:)));
K = 23;
F1121(idx(:,2,:)) = (1./(K - 3)).*((g(idx(:,2,:)).^2).*(K - 2).*F1123(idx(:,2,:)) + nu_eta(idx(:,2,:)).*E2121(idx(:,2,:)) - nu_xi(idx(:,2,:)).*E1221(idx(:,2,:)));
K = 21;
F1119(idx(:,2,:)) = (1./(K - 3)).*((g(idx(:,2,:)).^2).*(K - 2).*F1121(idx(:,2,:)) + nu_eta(idx(:,2,:)).*E2119(idx(:,2,:)) - nu_xi(idx(:,2,:)).*E1219(idx(:,2,:)));
K = 19;
F1117(idx(:,2,:)) = (1./(K - 3)).*((g(idx(:,2,:)).^2).*(K - 2).*F1119(idx(:,2,:)) + nu_eta(idx(:,2,:)).*E2117(idx(:,2,:)) - nu_xi(idx(:,2,:)).*E1217(idx(:,2,:)));
K = 17;
F1115(idx(:,2,:)) = (1./(K - 3)).*((g(idx(:,2,:)).^2).*(K - 2).*F1117(idx(:,2,:)) + nu_eta(idx(:,2,:)).*E2115(idx(:,2,:)) - nu_xi(idx(:,2,:)).*E1215(idx(:,2,:)));
K = 15;
F1113(idx(:,2,:)) = (1./(K - 3)).*((g(idx(:,2,:)).^2).*(K - 2).*F1115(idx(:,2,:)) + nu_eta(idx(:,2,:)).*E2113(idx(:,2,:)) - nu_xi(idx(:,2,:)).*E1213(idx(:,2,:)));
K = 13;
F1111(idx(:,2,:)) = (1./(K - 3)).*((g(idx(:,2,:)).^2).*(K - 2).*F1113(idx(:,2,:)) + nu_eta(idx(:,2,:)).*E2111(idx(:,2,:)) - nu_xi(idx(:,2,:)).*E1211(idx(:,2,:)));
K = 11;
F119(idx(:,2,:)) = (1./(K - 3)).*((g(idx(:,2,:)).^2).*(K - 2).*F1111(idx(:,2,:)) + nu_eta(idx(:,2,:)).*E219(idx(:,2,:)) - nu_xi(idx(:,2,:)).*E129(idx(:,2,:)));
K = 9;
F117(idx(:,2,:)) = (1./(K - 3)).*((g(idx(:,2,:)).^2).*(K - 2).*F119(idx(:,2,:)) + nu_eta(idx(:,2,:)).*E217(idx(:,2,:)) - nu_xi(idx(:,2,:)).*E127(idx(:,2,:)));
K = 7;
F115(idx(:,2,:)) = (1./(K - 3)).*((g(idx(:,2,:)).^2).*(K - 2).*F117(idx(:,2,:)) + nu_eta(idx(:,2,:)).*E215(idx(:,2,:)) - nu_xi(idx(:,2,:)).*E125(idx(:,2,:)));
K = 5;
F113(idx(:,2,:)) = (1./(K - 3)).*((g(idx(:,2,:)).^2).*(K - 2).*F115(idx(:,2,:)) + nu_eta(idx(:,2,:)).*E213(idx(:,2,:)) - nu_xi(idx(:,2,:)).*E123(idx(:,2,:)));

% PAGE 130
% 3.)
% a)
% (i)
idx = abs(nu_eta) <= abs(nu_xi);
F121 = zeros(size(F111));
F131 = zeros(size(F111));
F141 = zeros(size(F111));
F151 = zeros(size(F111));
F211 = zeros(size(F111));
F221 = zeros(size(F111));
F231 = zeros(size(F111));
F241 = zeros(size(F111));
F311 = zeros(size(F111));
F321 = zeros(size(F111));
F331 = zeros(size(F111));
F411 = zeros(size(F111));
F421 = zeros(size(F111));
F511 = zeros(size(F111));

h = repmat(h,1,1,3);
F121(idx) = (a(idx).*nu_eta(idx).*F111(idx) + nu_xi(idx).*E11n1(idx));
N = 3;
F131(idx) = (1./(N - 1)).*((2.*N - 3).*a(idx).*nu_eta(idx).*F121(idx) - (a(idx).^2 + (nu_xi(idx).^2).*(h(idx).^2)).*F111(idx) + nu_xi(idx).*E12n1(idx));
N = 4;
F141(idx) = (1./(N - 1)).*((2.*N - 3).*a(idx).*nu_eta(idx).*F131(idx) - (a(idx).^2 + (nu_xi(idx).^2).*(h(idx).^2)).*F121(idx) + nu_xi(idx).*E13n1(idx));
N = 5;
F151(idx) = (1./(N - 1)).*((2.*N - 3).*a(idx).*nu_eta(idx).*F141(idx) - (a(idx).^2 + (nu_xi(idx).^2).*(h(idx).^2)).*F131(idx) + nu_xi(idx).*E14n1(idx));

% PAGE 131
% (ii)
F211(idx) = (-nu_eta(idx)./nu_xi(idx)).*F121(idx) + (a(idx)./nu_xi(idx)).*F111(idx);
F221(idx) = (-nu_eta(idx)./nu_xi(idx)).*F131(idx) + (a(idx)./nu_xi(idx)).*F121(idx);
F231(idx) = (-nu_eta(idx)./nu_xi(idx)).*F141(idx) + (a(idx)./nu_xi(idx)).*F131(idx);
F241(idx) = (-nu_eta(idx)./nu_xi(idx)).*F151(idx) + (a(idx)./nu_xi(idx)).*F141(idx);

F311(idx) = (-nu_eta(idx)./nu_xi(idx)).*F221(idx) + (a(idx)./nu_xi(idx)).*F211(idx);
F321(idx) = (-nu_eta(idx)./nu_xi(idx)).*F231(idx) + (a(idx)./nu_xi(idx)).*F221(idx);
F331(idx) = (-nu_eta(idx)./nu_xi(idx)).*F241(idx) + (a(idx)./nu_xi(idx)).*F231(idx);

F411(idx) = (-nu_eta(idx)./nu_xi(idx)).*F321(idx) + (a(idx)./nu_xi(idx)).*F311(idx);
F421(idx) = (-nu_eta(idx)./nu_xi(idx)).*F331(idx) + (a(idx)./nu_xi(idx)).*F321(idx);

F511(idx) = (-nu_eta(idx)./nu_xi(idx)).*F421(idx) + (a(idx)./nu_xi(idx)).*F411(idx);

% b)
% (i)
idx = abs(nu_xi) <= abs(nu_eta);

M = 2;
F211(idx) = (1./(M - 1)).*((2.*M - 3).*a(idx).*nu_xi(idx).*F111(idx) - nu_eta(idx).*E11n1(idx));
M = 3;
F311(idx) = (1./(M - 1)).*((2.*M - 3).*a(idx).*nu_xi(idx).*F211(idx) - (M - 2).*(a(idx).^2 + (nu_eta(idx).^2).*(h(idx).^2)).*F111(idx) - nu_eta(idx).*E21n1(idx));
M = 4;
F411(idx) = (1./(M - 1)).*((2.*M - 3).*a(idx).*nu_xi(idx).*F311(idx) - (M - 2).*(a(idx).^2 + (nu_eta(idx).^2).*(h(idx).^2)).*F211(idx) - nu_eta(idx).*E31n1(idx));
M = 5;
F511(idx) = (1./(M - 1)).*((2.*M - 3).*a(idx).*nu_xi(idx).*F411(idx) - (M - 2).*(a(idx).^2 + (nu_eta(idx).^2).*(h(idx).^2)).*F311(idx) - nu_eta(idx).*E41n1(idx));

% (ii)
F121(idx) = (-nu_xi(idx)./nu_eta(idx)).*F211(idx) + (a(idx)./nu_eta(idx)).*F111(idx);
F221(idx) = (-nu_xi(idx)./nu_eta(idx)).*F311(idx) + (a(idx)./nu_eta(idx)).*F211(idx);
F321(idx) = (-nu_xi(idx)./nu_eta(idx)).*F411(idx) + (a(idx)./nu_eta(idx)).*F311(idx);
F421(idx) = (-nu_xi(idx)./nu_eta(idx)).*F511(idx) + (a(idx)./nu_eta(idx)).*F411(idx);

F131(idx) = (-nu_xi(idx)./nu_eta(idx)).*F221(idx) + (a(idx)./nu_eta(idx)).*F121(idx);
F231(idx) = (-nu_xi(idx)./nu_eta(idx)).*F321(idx) + (a(idx)./nu_eta(idx)).*F221(idx);
F331(idx) = (-nu_xi(idx)./nu_eta(idx)).*F421(idx) + (a(idx)./nu_eta(idx)).*F321(idx);

F141(idx) = (-nu_xi(idx)./nu_eta(idx)).*F231(idx) + (a(idx)./nu_eta(idx)).*F131(idx);
F241(idx) = (-nu_xi(idx)./nu_eta(idx)).*F331(idx) + (a(idx)./nu_eta(idx)).*F231(idx);

F151(idx) = (-nu_xi(idx)./nu_eta(idx)).*F241(idx) + (a(idx)./nu_eta(idx)).*F141(idx);

% PAGE 132
% 4.)
F123 = nu_eta.*a.*F113 - nu_xi.*E111;

% 5.)
F133 = 2.*a.*nu_eta.*F123 - (a.^2 + (nu_xi.^2).*(h.^2)).*F113 + (nu_xi.^2).*F111;
F143 = 2.*a.*nu_eta.*F133 - (a.^2 + (nu_xi.^2).*(h.^2)).*F123 + (nu_xi.^2).*F121;
F153 = 2.*a.*nu_eta.*F143 - (a.^2 + (nu_xi.^2).*(h.^2)).*F133 + (nu_xi.^2).*F131;

% Correct F for proximity to perimeter

%% H

% abs(h) >= delta_h
h = h(:,:,1);

% PAGE 124
% 1.)
H111 = -abs(h).*sum(atan2(a.*(L2.*c1 - L1.*c2), c1.*c2 + (a.^2).*L1.*L2),3) + sum(a.*F111, 3);

% PAGE 125
% 2.)
idx = [abs(h) >= del_h.*d_H, abs(h) < del_h.*d_H & ~idx_on_element abs(h) < del_h.*d_H & idx_on_element]; 
H113(idx(:,1),1) = (1./h(idx(:,1)).^2).*(-1.*H111(idx(:,1)) + sum(a(idx(:,1),:,:).*F111(idx(:,1),:,:), 3));
H115(idx(:,1),1) = (1./(3.*h(idx(:,1)).^2)).*(H113(idx(:,1)) + sum(a(idx(:,1),:,:).*F113(idx(:,1),:,:), 3));

% PAGE 126
% 1.)
H1121 = zeros(size(H111));
% 2.)
K = 21;
H1119(idx(:,2),1) = (1./(K - 4)).*((h(idx(:,2)).^2).*(K - 2).*H1121(idx(:,2)) - sum(a(idx(:,2),:,:).*F1119(idx(:,2),:,:), 3)); 
K = 19;
H1117(idx(:,2),1) = (1./(K - 4)).*((h(idx(:,2)).^2).*(K - 2).*H1119(idx(:,2)) - sum(a(idx(:,2),:,:).*F1117(idx(:,2),:,:), 3));
K = 17;
H1115(idx(:,2),1) = (1./(K - 4)).*((h(idx(:,2)).^2).*(K - 2).*H1117(idx(:,2)) - sum(a(idx(:,2),:,:).*F1115(idx(:,2),:,:), 3));
K = 15;
H1113(idx(:,2),1) = (1./(K - 4)).*((h(idx(:,2)).^2).*(K - 2).*H1115(idx(:,2)) - sum(a(idx(:,2),:,:).*F1113(idx(:,2),:,:), 3));
K = 13;
H1111(idx(:,2),1) = (1./(K - 4)).*((h(idx(:,2)).^2).*(K - 2).*H1113(idx(:,2)) - sum(a(idx(:,2),:,:).*F1111(idx(:,2),:,:), 3));
K = 11;
H119(idx(:,2),1) = (1./(K - 4)).*((h(idx(:,2)).^2).*(K - 2).*H1111(idx(:,2)) - sum(a(idx(:,2),:,:).*F119(idx(:,2),:,:), 3));
K = 9;
H117(idx(:,2),1) = (1./(K - 4)).*((h(idx(:,2)).^2).*(K - 2).*H119(idx(:,2)) - sum(a(idx(:,2),:,:).*F117(idx(:,2),:,:), 3));
K = 7;
H115(idx(:,2),1) = (1./(K - 4)).*((h(idx(:,2)).^2).*(K - 2).*H117(idx(:,2)) - sum(a(idx(:,2),:,:).*F115(idx(:,2),:,:), 3));
K = 5;
H113(idx(:,2),1) = (1./(K - 4)).*((h(idx(:,2)).^2).*(K - 2).*H115(idx(:,2)) - sum(a(idx(:,2),:,:).*F113(idx(:,2),:,:), 3));
K = 3;
H111(idx(:,2),1) = (1./(K - 4)).*((h(idx(:,2)).^2).*(K - 2).*H113(idx(:,2)) - sum(a(idx(:,2),:,:).*F111(idx(:,2),:,:), 3));

% PAGE 125
% 3.)
H211 = (1/2).*((h.^2).*sum(nu_xi.*F111,3) + sum(a.*F211,3));
H221 = (1/3).*((h.^2).*sum(nu_xi.*F121,3) + sum(a.*F221,3));
H231 = (1/3).*((h.^2).*sum(nu_xi.*F131,3) + sum(a.*F231,3));
H241 = (1/3).*((h.^2).*sum(nu_xi.*F141,3) + sum(a.*F241,3));

% 4.)
H121 = (1/2).*((h.^2).*sum(nu_eta.*F111,3) + sum(a.*F121,3));
H131 = (1/3).*(-(h.^2).*H111 + (h.^2).*sum(nu_eta.*F121,3) + sum(a.*F131,3));
H141 = (1/4).*(-(h.^2).*2.*H121 + (h.^2).*sum(nu_eta.*F131,3) + sum(a.*F141,3));
H151 = (1/5).*(-(h.^2).*2.*H131 + (h.^2).*sum(nu_eta.*F141,3) + sum(a.*F151,3));

% 5.)
M = 3; N = 1;
H311 = (1./(M + N - 1)).*(-(h.^2).*(M - 2).*H111 + (h.^2).*sum(nu_xi.*F211,3) + sum(a.*F311,3));
M = 4; N = 1;
H411 = (1./(M + N - 1)).*(-(h.^2).*(M - 2).*H211 + (h.^2).*sum(nu_xi.*F311,3) + sum(a.*F411,3));
M = 5; N = 1;
H511 = (1./(M + N - 1)).*(-(h.^2).*(M - 2).*H311 + (h.^2).*sum(nu_xi.*F411,3) + sum(a.*F511,3));

M = 3; N = 2;
H321 = (1./(M + N - 1)).*(-(h.^2).*(M - 2).*H121 + (h.^2).*sum(nu_xi.*F221,3) + sum(a.*F321,3));
M = 4; N = 2;
H421 = (1./(M + N - 1)).*(-(h.^2).*(M - 2).*H221 + (h.^2).*sum(nu_xi.*F321,3) + sum(a.*F421,3));

M = 3; N = 3;
H331 = (1./(M + N - 1)).*(-(h.^2).*(M - 2).*H131 + (h.^2).*sum(nu_xi.*F231,3) + sum(a.*F331,3));


% PAGE 126
% 6.)
K = 3;
H123 = (1./(K - 2)).*(-sum(nu_eta.*F111,3));
N = 3; K = 3;
H133 = (1./(K - 2)).*((N - 2).*H111 - sum(nu_eta.*F121,3));
N = 4; K = 3;
H143 = (1./(K - 2)).*((N - 2).*H121 - sum(nu_eta.*F131,3));
N = 5; K = 3;
H153 = (1./(K - 2)).*((N - 2).*H131 - sum(nu_eta.*F141,3));

K = 5;
H125 = (1./(K - 2)).*(-sum(nu_eta.*F113,3));
N = 3; K = 5;
H135 = (1./(K - 2)).*((N - 2).*H113 - sum(nu_eta.*F123,3));
N = 4; K = 5;
H145 = (1./(K - 2)).*((N - 2).*H123 - sum(nu_eta.*F133,3));
N = 5; K = 5;
H155 = (1./(K - 2)).*((N - 2).*H133 - sum(nu_eta.*F143,3));

% 7.)
K = 3;
H213 = (1./(K - 2)).*(-sum(nu_xi.*F111,3));
H223 = (1./(K - 2)).*(-sum(nu_xi.*F121,3));
H233 = (1./(K - 2)).*(-sum(nu_xi.*F131,3));
H243 = (1./(K - 2)).*(-sum(nu_xi.*F141,3));

H215 = (1./(K - 2)).*(-sum(nu_xi.*F113,3));
H225 = (1./(K - 2)).*(-sum(nu_xi.*F123,3));
H235 = (1./(K - 2)).*(-sum(nu_xi.*F133,3));
H245 = (1./(K - 2)).*(-sum(nu_xi.*F143,3));

% 8.)
H313 = -H133 - (h.^2).*H113 + H111;
H413 = -H233 - (h.^2).*H213 + H211;
H323 = -H143 - (h.^2).*H123 + H121;
H333 = -H153 - (h.^2).*H133 + H131;
H513 = -H333 - (h.^2).*H313 + H311;

H315 = -H135 - (h.^2).*H115 + H113;
H415 = -H235 - (h.^2).*H215 + H213;
H325 = -H145 - (h.^2).*H125 + H123;
H335 = -H155 - (h.^2).*H135 + H133;
H515 = -H335 - (h.^2).*H315 + H313;


%%
J11 = [reshape(3.*h.*H215,1,1,[]); reshape(3.*h.*H125,1,1,[]); reshape(H113 - 3.*(h.^2).*H115,1,1,[])];
J12 = [reshape(3.*h.*H225,1,1,[]); reshape(3.*h.*H135,1,1,[]); reshape(H123 - 3.*(h.^2).*H125,1,1,[])];
J13 = [reshape(3.*h.*H235,1,1,[]); reshape(3.*h.*H145,1,1,[]); reshape(H133 - 3.*(h.^2).*H135,1,1,[])];
J21 = [reshape(3.*h.*H315,1,1,[]); reshape(3.*h.*H225,1,1,[]); reshape(H213 - 3.*(h.^2).*H215,1,1,[])];
J22 = [reshape(3.*h.*H325,1,1,[]); reshape(3.*h.*H235,1,1,[]); reshape(H223 - 3.*(h.^2).*H225,1,1,[])];
J31 = [reshape(3.*h.*H415,1,1,[]); reshape(3.*h.*H325,1,1,[]); reshape(H313 - 3.*(h.^2).*H315,1,1,[])];

%%
y = reshape(y_m,1,1,[]);
x = reshape(x_m,1,1,[]);
infl_loc(:,1,:) = 0.5.*(y.^2).*J11 + y.*J12 + 0.5.*J13;
infl_loc(:,2,:) = J11.*y + J12;
infl_loc(:,3,:) = 0.5.*(x.^2).*J11 + x.*J21 + 0.5.*J31;
infl_loc(:,4,:) = J11.*x + J21;
infl_loc(:,5,:) = J11.*y.*x + J12.*x + J21.*y + J22;
infl_loc(:,6,:) = J11;

infl_loc(:,:,idx_on_element) = infl_loc(:,:,idx_on_element).*0;
disp('Turning off on element');

end