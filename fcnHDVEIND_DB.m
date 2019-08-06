function [infl_loc] = fcnHDVEIND_DB(dvenum, dvetype, fpg, matPLEX, matROTANG, matCONTROL, vecBI, ztol)
warning('on')

fpl = fcnGLOBSTAR(fpg - matCONTROL(dvenum,:), matROTANG(dvenum,:));
len = size(fpl,1);

x_m = fpl(:,1);
y_m = fpl(:,2);

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
idx_on_element = y_m >= te_eta - margin & y_m <= le_eta + margin & x_m >= xi_left - margin & x_m <= xi_right + margin & abs(fpl(:,3)) < ztol*2;

margin = 1e-5;
idx_flp = xi_3 < xi_1; % Flipping influence of elements that need a good flippin
if any(abs(xi_2 - xi_3) < margin & abs(xi_1 - xi_2) > margin)
    disp('Issue in element orientation in HDVEIND.');
end

%%
h = fpl(:,3);
h(idx_on_element) = sign(h(idx_on_element)).*sqrt(fpl(idx_on_element,3).^2 + ztol.^2);
h(idx_on_element & h == 0) = ztol;
abs_h = abs(h);
hs = h.^2;

k = zeros(size(h,1),1,3);
% k(idx_on_element,:,:) = ztol;
% k = repmat(ztol, size(h,1), 1, 3);

%%
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

% tmp = nan(len,1,3);
% idx = sign(L1) ~= sign(L2);
% tmp(idx) = a(idx);
% idx = (sign(L1) == sign(L2)) & ~idx;
% tmp(idx) = max([abs(a(idx)) min([abs(L1(idx)) abs(L2(idx))],[],2)], [], 2);
% if any(isnan(tmp),'all')
%     disp('Issue in determining d_H');
% end
% d_H = min(abs(tmp),[],3);

g = sqrt(a.^2 + hs);
c1 = g.^2 + abs_h.*sqrt(L1.^2 + g.^2);
c2 = g.^2 + abs_h.*sqrt(L2.^2 + g.^2);

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
idx = (L1 < 0 & L2 >= 0) | (L1 >= 0 & L2 < 0);
F111(idx) = log(((rho1(idx) - L1(idx)).*(rho2(idx) + L2(idx)) + k(idx))./(g(idx).^2 + k(idx)));

% idx = [g > del_h.*d_H, g <= del_h.*d_H];
idx = [g > 1e-5, g <= 1e-5];
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
F211 = zeros(size(F111));

F121(idx) = (a(idx).*nu_eta(idx).*F111(idx) + nu_xi(idx).*E11n1(idx));

% PAGE 131
% (ii)
F211(idx) = (-nu_eta(idx)./nu_xi(idx)).*F121(idx) + (a(idx)./nu_xi(idx)).*F111(idx);

% b)
% (i)
idx = abs(nu_xi) <= abs(nu_eta);

M = 2;
F211(idx) = (1./(M - 1)).*((2.*M - 3).*a(idx).*nu_xi(idx).*F111(idx) - nu_eta(idx).*E11n1(idx));

% (ii)
F121(idx) = (-nu_xi(idx)./nu_eta(idx)).*F211(idx) + (a(idx)./nu_eta(idx)).*F111(idx);

% PAGE 132
% 4.)
F123 = nu_eta.*a.*F113 - nu_xi.*E111;

% 5.)
F133 = 2.*a.*nu_eta.*F123 - (a.^2 + (nu_xi.^2).*(hs)).*F113 + (nu_xi.^2).*F111;

%% H

% PAGE 124
% 1.)
H111 = -abs_h.*sum(atan2(a.*(L2.*c1 - L1.*c2), c1.*c2 + (a.^2).*L1.*L2),3) + sum(a.*F111, 3);

% PAGE 125
% 2.)
% idx = [abs_h > del_h.*d_H, abs_h <= del_h.*d_H]; 
idx = [abs_h > 1e-5, abs_h <= 1e-5]; 

H113(idx(:,1),1) = (1./hs(idx(:,1))).*(-1.*H111(idx(:,1)) + sum(a(idx(:,1),:,:).*F111(idx(:,1),:,:), 3));
H115(idx(:,1),1) = (1./(3.*hs(idx(:,1)))).*(H113(idx(:,1)) + sum(a(idx(:,1),:,:).*F113(idx(:,1),:,:), 3));

% PAGE 126
% 1.)
H1121 = zeros(size(H111));

% 2.)
K = 21;
H1119(idx(:,2),1) = (1./(K - 4)).*(hs(idx(:,2)).*(K - 2).*H1121(idx(:,2)) - sum(a(idx(:,2),:,:).*F1119(idx(:,2),:,:), 3)); 
K = 19;
H1117(idx(:,2),1) = (1./(K - 4)).*(hs(idx(:,2)).*(K - 2).*H1119(idx(:,2)) - sum(a(idx(:,2),:,:).*F1117(idx(:,2),:,:), 3));
K = 17;
H1115(idx(:,2),1) = (1./(K - 4)).*(hs(idx(:,2)).*(K - 2).*H1117(idx(:,2)) - sum(a(idx(:,2),:,:).*F1115(idx(:,2),:,:), 3));
K = 15;
H1113(idx(:,2),1) = (1./(K - 4)).*(hs(idx(:,2)).*(K - 2).*H1115(idx(:,2)) - sum(a(idx(:,2),:,:).*F1113(idx(:,2),:,:), 3));
K = 13;
H1111(idx(:,2),1) = (1./(K - 4)).*(hs(idx(:,2)).*(K - 2).*H1113(idx(:,2)) - sum(a(idx(:,2),:,:).*F1111(idx(:,2),:,:), 3));
K = 11;
H119(idx(:,2),1) = (1./(K - 4)).*(hs(idx(:,2)).*(K - 2).*H1111(idx(:,2)) - sum(a(idx(:,2),:,:).*F119(idx(:,2),:,:), 3));
K = 9;
H117(idx(:,2),1) = (1./(K - 4)).*(hs(idx(:,2)).*(K - 2).*H119(idx(:,2)) - sum(a(idx(:,2),:,:).*F117(idx(:,2),:,:), 3));
K = 7;
H115(idx(:,2),1) = (1./(K - 4)).*(hs(idx(:,2)).*(K - 2).*H117(idx(:,2)) - sum(a(idx(:,2),:,:).*F115(idx(:,2),:,:), 3));
K = 5;
H113(idx(:,2),1) = (1./(K - 4)).*(hs(idx(:,2)).*(K - 2).*H115(idx(:,2)) - sum(a(idx(:,2),:,:).*F113(idx(:,2),:,:), 3));
K = 3;
H111(idx(:,2),1) = (1./(K - 4)).*(hs(idx(:,2)).*(K - 2).*H113(idx(:,2)) - sum(a(idx(:,2),:,:).*F111(idx(:,2),:,:), 3));

% PAGE 127
% idx = idx(:,2) & idx_on_element & abs(h) < (1e-8)*d_H;
% idx = []
% H111(idx) = H111(idx) - eps(1,1,1,h(idx));
% H113(idx) = H113(idx) - eps(1,1,3,h(idx));
% H115(idx) = H115(idx) - eps(1,1,5,h(idx));

% PAGE 126
% 6.)
K = 3;
H123 = (1./(K - 2)).*(-sum(nu_eta.*F111,3));
N = 3; K = 3;
H133 = (1./(K - 2)).*((N - 2).*H111 - sum(nu_eta.*F121,3));

K = 5;
H125 = (1./(K - 2)).*(-sum(nu_eta.*F113,3));
N = 3;
H135 = (1./(K - 2)).*((N - 2).*H113 - sum(nu_eta.*F123,3));
N = 4;
H145 = (1./(K - 2)).*((N - 2).*H123 - sum(nu_eta.*F133,3));

% 7.)
K = 3;
H213 = (1./(K - 2)).*(-sum(nu_xi.*F111,3));
H223 = (1./(K - 2)).*(-sum(nu_xi.*F121,3));

K = 5;
H215 = (1./(K - 2)).*(-sum(nu_xi.*F113,3));
H225 = (1./(K - 2)).*(-sum(nu_xi.*F123,3));
H235 = (1./(K - 2)).*(-sum(nu_xi.*F133,3));

% 8.)
H313 = -H133 - (hs).*H113 + H111;
H315 = -H135 - (hs).*H115 + H113;
H415 = -H235 - (hs).*H215 + H213;
H325 = -H145 - (hs).*H125 + H123;

%%
h = fpl(:,3);
J11 = [reshape(3.*h.*H215,1,1,[]); reshape(3.*h.*H125,1,1,[]); reshape(H113 - 3.*(hs).*H115,1,1,[])];
J12 = [reshape(3.*h.*H225,1,1,[]); reshape(3.*h.*H135,1,1,[]); reshape(H123 - 3.*(hs).*H125,1,1,[])];
J13 = [reshape(3.*h.*H235,1,1,[]); reshape(3.*h.*H145,1,1,[]); reshape(H133 - 3.*(hs).*H135,1,1,[])];
J21 = [reshape(3.*h.*H315,1,1,[]); reshape(3.*h.*H225,1,1,[]); reshape(H213 - 3.*(hs).*H215,1,1,[])];
J22 = [reshape(3.*h.*H325,1,1,[]); reshape(3.*h.*H235,1,1,[]); reshape(H223 - 3.*(hs).*H225,1,1,[])];
J31 = [reshape(3.*h.*H415,1,1,[]); reshape(3.*h.*H325,1,1,[]); reshape(H313 - 3.*(hs).*H315,1,1,[])];

%%
x = reshape(x_m,1,1,[]);
y = reshape(y_m,1,1,[]);

infl_loc = nan(3,6,length(dvenum));
infl_loc(:,1,:) = 0.5.*(y.^2).*J11 + y.*J12 + 0.5.*J13;
infl_loc(:,2,:) = J11.*y + J12;
infl_loc(:,3,:) = 0.5.*(x.^2).*J11 + x.*J21 + 0.5.*J31;
infl_loc(:,4,:) = J11.*x + J21;
infl_loc(:,5,:) = J11.*y.*x + J12.*x + J21.*y + J22;
infl_loc(:,6,:) = J11;

if any(vecBI)
    infl_BI = fcnBOUNDIND(dvenum(vecBI), dvetype(vecBI), fpg(vecBI,:), matPLEX, matROTANG, matCONTROL, vecBI(vecBI), 0);
    infl_loc(:,:,vecBI) = infl_loc(:,:,vecBI) - infl_BI;
end

infl_loc(:,:,idx_flp) = -infl_loc(:,:,idx_flp);

end




