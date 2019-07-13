% clc
clear

load('single_dve.mat');

% a -> (or a_bar) distance from edge, shortest distance, perpendicular to edge
% h -> height above surface, z_m
% L -> distance (from pt1 or pt2) along edge to field point intersection

% fpg = [1 2.2 1; 0 0 0; 0 1 0; -1 -1 -0.4];
fpg = [1.3 2.2 2];
dvenum = ones(size(fpg,1),1);

%%
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

hFig1 = figure(1);
clf(1);
patch([xi_1(1) xi_2(1) xi_3(1)], [eta_1(1) eta_2(1) eta_3(1)],'w');
grid minor
box on
axis tight
axis equal
xlabel('Xi-Dir');
ylabel('Eta-Dir');
zlabel('Zi-Dir');
for i = 1:size(fpl,1)
    hold on
    scatter3(fpl(i,1), fpl(i,2), fpl(i,3), 'rs')
    hold off
end

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

g = sqrt(a.^2 + h.^2);
c1 = g.^2 + abs(h).*sqrt(L1.^2 + g.^2);
c2 = g.^2 + abs(h).*sqrt(L2.^2 + g.^2);

%% E
% MXK = 5, MXQ = 5, MXFK = 3;
% Out of plane
% E213, E215
% E123, E125

rho1 = sqrt(L1.^2 + g.^2);
rho2 = sqrt(L2.^2 + g.^2);

% a.)
E211 = (epts(:,1,:,2)./rho2) - (epts(:,1,:,1)./rho1); % I = 1

% b.)
E121 = (epts(:,2,:,2)./rho2) - (epts(:,2,:,1)./rho1); % I = 1

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

% Additional In plane (MXFK = 16 + 5 - 2 = 19):
% E217, E219, E2111, E2113, E2115, E2117, E2119
% E127, E129, E1211, E1213, E1215, E1217, E1219

%% F
% MXK = 5, MXQ = 5, MXFK = 3;
% Out of plane:
% F111, F113, F151, F241, F331, F421, F511, F123, F133, F143, F153

% Additional In plane (MXFK = 16 + 5 - 2 = 19):
% F151, F171, F191, F1111, F1131, F1151, F1171, F1191

% 1.)
F111 = zeros(size(a));

idx = L1 >= 0 & L2 >= 0;
F111(idx) = log((rho2(idx) + L2(idx))./(rho1(idx) + L1(idx)));

idx = L1 < 0 & L2 < 0;
F111(idx) = log((rho1(idx) - L1(idx))./(rho2(idx) - L2(idx)));

idx = L1 < 0 & L2 >= 0;
F111(idx) = log(((rho1(idx) - L1(idx)).*(rho2(idx) + L2(idx)))./(g(idx).^2));

% 2.)
F113 = (1./(g.^2)).*(-nu_n(:,2,:).*E211 + nu_n(:,1,:).*E121);

% 3.)
nu_xi = nu_n(:,1,:);
nu_eta = nu_n(:,2,:);

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
H113 = (1./h.^2).*(-1.*H111 + sum(a.*F111, 3));
H115 = (1./(3.*h.^2)).*(H113 + sum(a.*F113, 3));

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









