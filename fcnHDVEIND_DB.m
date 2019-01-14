function [infl_loc] = fcnHDVEIND_DB(dvenum, dvetype, fpg, matPLEX, matROTANG, matCONTROL)
warning('on')
% tol = 1e-10;
% tol2 = 1e-6;
% tol3 = 1e-5;

tol = 1e-10;
tol2 = 1e-4;
tol3 = 1e-4;


fpl = fcnGLOBSTAR(fpg - matCONTROL(dvenum,:), matROTANG(dvenum,:));

len = size(fpl,1);

x_m = fpl(:,1);
y_m = fpl(:,2);
z_m = fpl(:,3);

%% Checking state of field point with relation to element surface
margin_edge = 1e-5;
margin_on_element = 1e-5;

xi_1 = permute(matPLEX(1,1,dvenum),[3 2 1]);
xi_2 = permute(matPLEX(2,1,dvenum),[3 2 1]);
xi_3 = permute(matPLEX(3,1,dvenum),[3 2 1]);

eta_1 = permute(matPLEX(1,2,dvenum),[3 2 1]);
eta_2 = permute(matPLEX(2,2,dvenum),[3 2 1]);
eta_3 = permute(matPLEX(3,2,dvenum),[3 2 1]);

idx_flp = abs(xi_2 - xi_3) < margin_edge;
xi_tmp(idx_flp) = xi_3(idx_flp);
xi_3(idx_flp) = xi_1(idx_flp);
xi_1(idx_flp) = xi_tmp(idx_flp);
eta_tmp(idx_flp) = eta_3(idx_flp);
eta_3(idx_flp) = eta_1(idx_flp);
eta_1(idx_flp) = eta_tmp(idx_flp);

idx_rrg = eta_2 < eta_1;
eta_tmp(idx_rrg) = eta_2(idx_rrg);
eta_2(idx_rrg) = eta_1(idx_rrg);
eta_1(idx_rrg) = eta_tmp(idx_rrg);

% Checking which elements are on the element
C = (eta_3 - eta_2)./(xi_3 - xi_2);
D_LE = eta_2 - ((xi_2.*(eta_3 - eta_2))./(xi_3 - xi_2));
E = (eta_3 - eta_1)./(xi_3 - xi_1);
D_TE = eta_1 - ((xi_1.*(eta_3 - eta_1))./(xi_3 - xi_1));

le_eta = C.*x_m + D_LE;
te_eta = E.*x_m + D_TE;

xi_left = min([xi_1, xi_3],[],2);
xi_right = max([xi_1, xi_3],[],2);
idx_in_plane = abs(z_m) <= tol2;
idx_on_element = y_m >= te_eta - margin_edge & y_m <= le_eta + margin_edge & x_m >= xi_left - margin_edge & x_m <= xi_right + margin_edge & abs(z_m) <= margin_on_element;
idx_on_edge =   (abs(y_m - te_eta) < margin_edge & (xi_left - margin_edge <= x_m & x_m <= xi_right + margin_edge) & abs(z_m) <= margin_on_element) | ...
    (abs(y_m - le_eta) < margin_edge & (xi_left - margin_edge <= x_m & x_m <= xi_right + margin_edge) & abs(z_m) <= margin_on_element) | ...
    (abs(x_m - xi_left) < margin_edge & (te_eta - margin_edge <= y_m & y_m <= le_eta + margin_edge) & abs(z_m) <= margin_on_element) | ...
    (abs(x_m - xi_right) < margin_edge & (te_eta - margin_edge <= y_m & y_m <= le_eta + margin_edge) & abs(z_m) <= margin_on_element);
% disp(['Edge calls: ', num2str(sum(idx_on_edge))]);
if any(idx_on_edge)
    tmp = [x_m(idx_on_edge) y_m(idx_on_edge) z_m(idx_on_edge)] + (1e-2.*(-[x_m(idx_on_edge) y_m(idx_on_edge) z_m(idx_on_edge)]./(sqrt(sum([x_m(idx_on_edge) y_m(idx_on_edge) z_m(idx_on_edge)].^2,2)))));
    x_m(idx_on_edge) = tmp(:,1);
    y_m(idx_on_edge) = tmp(:,2);
    z_m(idx_on_edge) = tmp(:,3);
end

%%
alpha = z_m.^2;
N = [-C.*x_m - D_LE + y_m, -E.*x_m - D_TE + y_m];
S = [C.^2 + 1, E.^2 + 1];
T = [2.*C.*N(:,1), 2.*E.*N(:,2)];
u = [N(:,1).^2 + alpha, N(:,2).^2 + alpha];
L = [C E];
F_lim(:,:,1) = x_m - xi_3;
F_lim(:,:,2) = x_m - xi_1;

tol_F = 1e-5;
tol_S = 1e-14;
tol_T = 1e-14;
tol_u = 1e-5;
tol_alpha = 1e-14;

% Correcting zero F limits if we are in the plane of the element
idx = abs(F_lim(:,:,1)) < tol_F & alpha(:,1) <= tol_alpha;
sgn = sign(F_lim(idx,:,1));
sgn(sign(F_lim(idx,:,1)) == 0) = sign(F_lim(sign(F_lim(idx,:,1)) == 0,:,2));
F_lim(idx,:,1) = sgn.*tol_F;

idx = abs(F_lim(:,:,2)) < tol_F & alpha <= tol_alpha;
sgn = sign(F_lim(idx,:,2));
sgn(sign(F_lim(idx,:,2)) == 0) = sign(F_lim(sign(F_lim(idx,:,2)) == 0,:,1));
F_lim(idx,:,2) = sgn.*tol_F;

F_lim = repmat(F_lim,1,2,1);
alpha = repmat(alpha,1,2,1);

% Identifying which cases are special
idx = [];
idx(:,:,1)  = abs(S - 1) >  tol_S & abs(T) >  tol_T & abs(u) >  tol_u &  alpha >  tol_alpha;
idx(:,:,2)  = abs(S - 1) >  tol_S & abs(T) >  tol_T & abs(u) >  tol_u &  alpha <= tol_alpha;
idx(:,:,3)  = abs(S - 1) >  tol_S & abs(T) >  tol_T & abs(u) <= tol_u &  alpha >  tol_alpha;
idx(:,:,4)  = abs(S - 1) >  tol_S & abs(T) >  tol_T & abs(u) <= tol_u &  alpha <= tol_alpha;
idx(:,:,5)  = abs(S - 1) >  tol_S & abs(T) <= tol_T & abs(u) >  tol_u &  alpha >  tol_alpha;
idx(:,:,6)  = abs(S - 1) >  tol_S & abs(T) <= tol_T & abs(u) >  tol_u &  alpha <= tol_alpha;
idx(:,:,7)  = abs(S - 1) >  tol_S & abs(T) <= tol_T & abs(u) <= tol_u &  alpha >  tol_alpha;
idx(:,:,8)  = abs(S - 1) >  tol_S & abs(T) <= tol_T & abs(u) <= tol_u &  alpha <= tol_alpha;
idx(:,:,9)  = abs(S - 1) <= tol_S & abs(T) >  tol_T & abs(u) >  tol_u &  alpha >  tol_alpha;
idx(:,:,10) = abs(S - 1) <= tol_S & abs(T) >  tol_T & abs(u) >  tol_u &  alpha <= tol_alpha;
idx(:,:,11) = abs(S - 1) <= tol_S & abs(T) >  tol_T & abs(u) <= tol_u &  alpha >  tol_alpha;
idx(:,:,12) = abs(S - 1) <= tol_S & abs(T) >  tol_T & abs(u) <= tol_u &  alpha <= tol_alpha;
idx(:,:,13) = abs(S - 1) <= tol_S & abs(T) <= tol_T & abs(u) >  tol_u &  alpha >  tol_alpha;
idx(:,:,14) = abs(S - 1) <= tol_S & abs(T) <= tol_T & abs(u) >  tol_u &  alpha <= tol_alpha;
idx(:,:,15) = abs(S - 1) <= tol_S & abs(T) <= tol_T & abs(u) <= tol_u &  alpha >  tol_alpha;
idx(:,:,16) = abs(S - 1) <= tol_S & abs(T) <= tol_T & abs(u) <= tol_u &  alpha <= tol_alpha;
idx(:,:,17) = abs(S - 1) >  tol_S & abs(T) >  tol_T & abs(u) >  tol_u & abs(u - alpha) <=  tol_u;
idx(:,:,18) = abs(S - 1) >  tol_S & abs(T) <= tol_T & abs(u) >  tol_u & abs(u - alpha) <=  tol_u;
idx(:,:,19) = abs(S - 1) <= tol_S & abs(T) >  tol_T & abs(u) >  tol_u & abs(u - alpha) <=  tol_u;
idx(:,:,20) = abs(S - 1) <= tol_S & abs(T) <= tol_T & abs(u) >  tol_u & abs(u - alpha) <=  tol_u;
% idx(:,:,17:20) = false(size(idx(:,:,1:4)));


K0 = fcnK0(S, T, u, alpha, F_lim, tol, idx);
K1 = fcnK1(S, T, u, alpha, F_lim, tol, idx);
K2 = fcnK2(S, T, u, alpha, F_lim, tol, idx);
K3 = fcnK3(S, T, u, alpha, F_lim, tol, idx);
K4 = fcnK4(S, T, u, alpha, F_lim, tol, idx);
K5 = fcnK5(S, T, u, alpha, F_lim, tol, idx);
K6 = fcnK6(S, T, u, alpha, F_lim, tol, idx);
K7 = fcnK7(S, T, u, alpha, F_lim, tol, idx);

F(:,1) = sum(([0.1e1./0.3e1,-0.1e1./0.3e1]).*([(2.*L.*S+L).*K3+(6.*N.*S-3.*N).*K2+(-3.*L.*alpha+6.*L.*u).*K1 + ( +N.*alpha+2.*u.*N ).*K0]), 2);
F(:,2) = sum(([-0.1e1./0.3e1,0.1e1./0.3e1]).*([(2.*L.*S+L).*K4+(6.*N.*S-3.*N-(2.*L.*S+L).*x_m).*K3+(-3.*L.*alpha+6.*L.*u-(6.*N.*S-3.*N).*x_m).*K2+(N.*alpha+2.*u.*N-(-3.*L.*alpha+6.*L.*u).*x_m).*K1 + ( -(N.*alpha+2.*u.*N).*x_m ).*K0]), 2);
F(:,3) = sum(([0.1e1./0.3e1,-0.1e1./0.3e1]).*([(2.*L.*S+L).*K5+(6.*N.*S-3.*N-2.*(2.*L.*S+L).*x_m).*K4+(-3.*L.*alpha+6.*L.*u-2.*(6.*N.*S-3.*N).*x_m+(2.*L.*S+L).*x_m.^2).*K3+(N.*alpha+2.*u.*N-2.*(-3.*L.*alpha+6.*L.*u).*x_m+(6.*N.*S-3.*N).*x_m.^2).*K2+(-2.*(N.*alpha+2.*u.*N).*x_m+(-3.*L.*alpha+6.*L.*u).*x_m.^2).*K1 + ( +(N.*alpha+2.*u.*N).*x_m.^2 ).*K0]), 2);
F(:,4) = sum(([0.1e1./0.3e1,-0.1e1./0.3e1]).*([K4+(2.*L.*S.*y_m+y_m.*L).*K3+(6.*N.*y_m.*S-3.*N.*y_m+2.*alpha).*K2+(-3.*L.*alpha.*y_m+6.*L.*u.*y_m).*K1 + ( +N.*y_m.*alpha+2.*y_m.*u.*N+alpha.^2 ).*K0]), 2);
F(:,5) = sum(([0.1e1./0.3e1,-0.1e1./0.3e1]).*([(L.*S-L).*K5+(3.*N.*S-3.*N+2.*y_m).*K4+(2.*L.*S.*y_m.^2+L.*S.*alpha+L.*y_m.^2-4.*L.*alpha+3.*L.*u).*K3+(6.*N.*S.*y_m.^2+3.*N.*S.*alpha-3.*N.*y_m.^2-4.*N.*alpha+N.*u+4.*alpha.*y_m).*K2+(-3.*L.*alpha.*y_m.^2+6.*L.*u.*y_m.^2-3.*L.*alpha.^2+3.*L.*alpha.*u).*K1 + ( +N.*alpha.*y_m.^2+2.*N.*u.*y_m.^2-N.*alpha.^2+N.*alpha.*u+2.*y_m.*alpha.^2 ).*K0]), 2);
F(:,6) = sum(([-0.1e1./0.3e1,0.1e1./0.3e1]).*([K5+(2.*L.*S.*y_m+y_m.*L-x_m).*K4+(6.*N.*y_m.*S-3.*N.*y_m+2.*alpha-(2.*L.*S.*y_m+y_m.*L).*x_m).*K3+(-3.*L.*alpha.*y_m+6.*L.*u.*y_m-(6.*N.*y_m.*S-3.*N.*y_m+2.*alpha).*x_m).*K2+(N.*y_m.*alpha+2.*y_m.*u.*N+alpha.^2-(-3.*L.*alpha.*y_m+6.*L.*u.*y_m).*x_m).*K1 + ( -(N.*y_m.*alpha+2.*y_m.*u.*N+alpha.^2).*x_m ).*K0]), 2);
F(:,7) = sum(([-0.1e1./0.3e1,0.1e1./0.3e1]).*([(2.*L.*S+L).*K6+(-3.*x_m.*(2.*L.*S+L)+6.*N.*S-3.*N).*K5+(3.*x_m.^2.*(2.*L.*S+L)-3.*x_m.*(6.*N.*S-3.*N)-3.*L.*alpha+6.*L.*u).*K4+(-x_m.^3.*(2.*L.*S+L)+3.*x_m.^2.*(6.*N.*S-3.*N)-3.*x_m.*(-3.*L.*alpha+6.*L.*u)+N.*alpha+2.*u.*N).*K3+(-x_m.^3.*(6.*N.*S-3.*N)+3.*x_m.^2.*(-3.*L.*alpha+6.*L.*u)-3.*x_m.*(N.*alpha+2.*u.*N)).*K2+(-x_m.^3.*(-3.*L.*alpha+6.*L.*u)+3.*x_m.^2.*(N.*alpha+2.*u.*N)).*K1 + ( -x_m.^3.*(N.*alpha+2.*u.*N) ).*K0]), 2);
F(:,8) = sum(([0.1e1./0.3e1,-0.1e1./0.3e1]).*([(3.*S-1).*K6+(3.*L.*S.*y_m+6.*L.*N-3.*L.*y_m).*K5+(9.*N.*S.*y_m-9.*N.*y_m+6.*S.*alpha+3.*y_m.^2-3.*alpha+3.*u).*K4+(2.*L.*S.*y_m.^3+3.*L.*S.*alpha.*y_m+L.*y_m.^3+12.*L.*N.*alpha-12.*L.*alpha.*y_m+9.*L.*u.*y_m).*K3+(6.*N.*S.*y_m.^3+9.*N.*S.*alpha.*y_m-3.*N.*y_m.^3-12.*N.*alpha.*y_m+3.*N.*u.*y_m+3.*S.*alpha.^2+6.*alpha.*y_m.^2-3.*alpha.^2+6.*alpha.*u).*K2+(-3.*L.*alpha.*y_m.^3+6.*L.*u.*y_m.^3+6.*L.*N.*alpha.^2-9.*L.*alpha.^2.*y_m+9.*L.*alpha.*u.*y_m).*K1 + ( +N.*alpha.*y_m.^3+2.*N.*u.*y_m.^3-3.*N.*alpha.^2.*y_m+3.*N.*alpha.*u.*y_m+3.*alpha.^2.*y_m.^2-alpha.^3+3.*alpha.^2.*u ).*K0]), 2);
F(:,9) = sum(([0.1e1./0.3e1,-0.1e1./0.3e1]).*([K6+(2.*L.*S.*y_m+y_m.*L-2.*x_m).*K5+(x_m.^2-2.*x_m.*(2.*L.*S.*y_m+y_m.*L)+6.*N.*y_m.*S-3.*N.*y_m+2.*alpha).*K4+(x_m.^2.*(2.*L.*S.*y_m+y_m.*L)-2.*x_m.*(6.*N.*y_m.*S-3.*N.*y_m+2.*alpha)-3.*L.*alpha.*y_m+6.*L.*u.*y_m).*K3+(x_m.^2.*(6.*N.*y_m.*S-3.*N.*y_m+2.*alpha)-2.*x_m.*(-3.*L.*alpha.*y_m+6.*L.*u.*y_m)+N.*y_m.*alpha+2.*y_m.*u.*N+alpha.^2).*K2+(x_m.^2.*(-3.*L.*alpha.*y_m+6.*L.*u.*y_m)-2.*x_m.*(N.*y_m.*alpha+2.*y_m.*u.*N+alpha.^2)).*K1 + ( +x_m.^2.*(N.*y_m.*alpha+2.*y_m.*u.*N+alpha.^2) ).*K0]), 2);
F(:,10) = sum(([-0.1e1./0.3e1,0.1e1./0.3e1]).*([(L.*S-L).*K6+(3.*N.*S-3.*N+2.*y_m-(L.*S-L).*x_m).*K5+(2.*L.*S.*y_m.^2+L.*S.*alpha+L.*y_m.^2-4.*L.*alpha+3.*L.*u-(3.*N.*S-3.*N+2.*y_m).*x_m).*K4+(6.*N.*S.*y_m.^2+3.*N.*S.*alpha-3.*N.*y_m.^2-4.*N.*alpha+N.*u+4.*alpha.*y_m-(2.*L.*S.*y_m.^2+L.*S.*alpha+L.*y_m.^2-4.*L.*alpha+3.*L.*u).*x_m).*K3+(-3.*L.*alpha.*y_m.^2+6.*L.*u.*y_m.^2-3.*L.*alpha.^2+3.*L.*alpha.*u-(6.*N.*S.*y_m.^2+3.*N.*S.*alpha-3.*N.*y_m.^2-4.*N.*alpha+N.*u+4.*alpha.*y_m).*x_m).*K2+(N.*alpha.*y_m.^2+2.*N.*u.*y_m.^2-N.*alpha.^2+N.*alpha.*u+2.*y_m.*alpha.^2-(-3.*L.*alpha.*y_m.^2+6.*L.*u.*y_m.^2-3.*L.*alpha.^2+3.*L.*alpha.*u).*x_m).*K1 + ( -(N.*alpha.*y_m.^2+2.*N.*u.*y_m.^2-N.*alpha.^2+N.*alpha.*u+2.*y_m.*alpha.^2).*x_m ).*K0]), 2);
F(:,11) = sum(([0.1e1./0.3e1,-0.1e1./0.3e1]).*([(2.*L.*S+L).*K7+(6.*N.*S-3.*N-4.*(2.*L.*S+L).*x_m).*K6+(-3.*L.*alpha+6.*L.*u-4.*(6.*N.*S-3.*N).*x_m+6.*(2.*L.*S+L).*x_m.^2).*K5+(N.*alpha+2.*u.*N-4.*(-3.*L.*alpha+6.*L.*u).*x_m+6.*(6.*N.*S-3.*N).*x_m.^2-4.*(2.*L.*S+L).*x_m.^3).*K4+(-4.*(N.*alpha+2.*u.*N).*x_m+6.*(-3.*L.*alpha+6.*L.*u).*x_m.^2-4.*(6.*N.*S-3.*N).*x_m.^3+(2.*L.*S+L).*x_m.^4).*K3+(6.*(N.*alpha+2.*u.*N).*x_m.^2-4.*(-3.*L.*alpha+6.*L.*u).*x_m.^3+(6.*N.*S-3.*N).*x_m.^4).*K2+(-4.*(N.*alpha+2.*u.*N).*x_m.^3+(-3.*L.*alpha+6.*L.*u).*x_m.^4).*K1 + ( +(N.*alpha+2.*u.*N).*x_m.^4 ).*K0]), 2);
F(:,12) = sum(([0.1e1./0.3e1,-0.1e1./0.3e1]).*([(L.*S-L).*K7+(-2.*x_m.*(L.*S-L)+3.*N.*S-3.*N+2.*y_m).*K6+(x_m.^2.*(L.*S-L)-2.*x_m.*(3.*N.*S-3.*N+2.*y_m)+2.*L.*S.*y_m.^2+L.*S.*alpha+L.*y_m.^2-4.*L.*alpha+3.*L.*u).*K5+(x_m.^2.*(3.*N.*S-3.*N+2.*y_m)-2.*x_m.*(2.*L.*S.*y_m.^2+L.*S.*alpha+L.*y_m.^2-4.*L.*alpha+3.*L.*u)+6.*N.*S.*y_m.^2+3.*N.*S.*alpha-3.*N.*y_m.^2-4.*N.*alpha+N.*u+4.*alpha.*y_m).*K4+(x_m.^2.*(2.*L.*S.*y_m.^2+L.*S.*alpha+L.*y_m.^2-4.*L.*alpha+3.*L.*u)-2.*x_m.*(6.*N.*S.*y_m.^2+3.*N.*S.*alpha-3.*N.*y_m.^2-4.*N.*alpha+N.*u+4.*alpha.*y_m)-3.*L.*alpha.*y_m.^2+6.*L.*u.*y_m.^2-3.*L.*alpha.^2+3.*L.*alpha.*u).*K3+(x_m.^2.*(6.*N.*S.*y_m.^2+3.*N.*S.*alpha-3.*N.*y_m.^2-4.*N.*alpha+N.*u+4.*alpha.*y_m)-2.*x_m.*(-3.*L.*alpha.*y_m.^2+6.*L.*u.*y_m.^2-3.*L.*alpha.^2+3.*L.*alpha.*u)+N.*alpha.*y_m.^2+2.*N.*u.*y_m.^2-N.*alpha.^2+N.*alpha.*u+2.*y_m.*alpha.^2).*K2+(x_m.^2.*(-3.*L.*alpha.*y_m.^2+6.*L.*u.*y_m.^2-3.*L.*alpha.^2+3.*L.*alpha.*u)-2.*x_m.*(N.*alpha.*y_m.^2+2.*N.*u.*y_m.^2-N.*alpha.^2+N.*alpha.*u+2.*y_m.*alpha.^2)).*K1 + ( +x_m.^2.*(N.*alpha.*y_m.^2+2.*N.*u.*y_m.^2-N.*alpha.^2+N.*alpha.*u+2.*y_m.*alpha.^2) ).*K0]), 2);
F(:,13) = sum(([-0.1e1./0.3e1,0.1e1./0.3e1]).*([(4.*L.*S-L).*K7+(12.*N.*S-12.*S.*y_m-9.*N+4.*y_m).*K6+(-6.*L.*S.*y_m.^2-24.*L.*N.*y_m+8.*L.*S.*alpha+6.*L.*y_m.^2-11.*L.*alpha+12.*L.*u).*K5+(-18.*N.*S.*y_m.^2+24.*N.*S.*alpha+18.*N.*y_m.^2-24.*S.*alpha.*y_m-4.*y_m.^3-19.*N.*alpha+4.*N.*u+12.*alpha.*y_m-12.*u.*y_m).*K4+(-2.*L.*S.*y_m.^4-6.*L.*S.*alpha.*y_m.^2-L.*y_m.^4-48.*L.*N.*alpha.*y_m+4.*L.*S.*alpha.^2+24.*L.*alpha.*y_m.^2-18.*L.*u.*y_m.^2-19.*L.*alpha.^2+24.*L.*alpha.*u).*K3+(-6.*N.*S.*y_m.^4-18.*N.*S.*alpha.*y_m.^2+3.*N.*y_m.^4+12.*N.*S.*alpha.^2+24.*N.*alpha.*y_m.^2-6.*N.*u.*y_m.^2-12.*S.*alpha.^2.*y_m-8.*alpha.*y_m.^3-11.*N.*alpha.^2+8.*N.*alpha.*u+12.*alpha.^2.*y_m-24.*alpha.*u.*y_m).*K2+(3.*L.*alpha.*y_m.^4-6.*L.*u.*y_m.^4-24.*L.*N.*alpha.^2.*y_m+18.*L.*alpha.^2.*y_m.^2-18.*L.*alpha.*u.*y_m.^2-9.*L.*alpha.^3+12.*L.*alpha.^2.*u).*K1 + ( -N.*alpha.*y_m.^4-2.*N.*u.*y_m.^4+6.*N.*alpha.^2.*y_m.^2-6.*N.*alpha.*u.*y_m.^2-4.*alpha.^2.*y_m.^3-N.*alpha.^3+4.*N.*alpha.^2.*u+4.*y_m.*alpha.^3-12.*alpha.^2.*u.*y_m ).*K0]), 2) ...
                - fcnH_7(E, -E.*x_m - D_TE + y_m, alpha(:,1), F_lim(:,1,1), F_lim(:,1,2), tol) + fcnH_7(C, -C.*x_m - D_LE + y_m, alpha(:,1), F_lim(:,1,1), F_lim(:,1,2), tol);
F = real(F);

%%
tmp11 = (-3/2).*F(:,5).*x_m.*z_m + (3/2).*F(:,10).*z_m;
tmp12 = -3.*F(:,4).*x_m.*z_m + 3.*F(:,6).*z_m;
tmp13 = (-3/2).*F(:,3).*x_m.*z_m + (3/2).*F(:,7).*z_m;
tmp14 = -3.*F(:,2).*x_m.*z_m + 3.*F(:,3).*z_m;
tmp15 = -3.*F(:,1).*x_m.*z_m + 3.*F(:,2).*z_m;

tmp21 = (3/2).*F(:,8).*z_m - (3/2).*F(:,5).*y_m.*z_m;
tmp22 = -3.*F(:,4).*y_m.*z_m + 3.*F(:,5).*z_m;
tmp23 = (3/2).*F(:,9).*z_m - (3/2).*F(:,3).*y_m.*z_m;
tmp24 = -3.*F(:,2).*y_m.*z_m + 3.*F(:,6).*z_m;
tmp25 = -3.*F(:,1).*y_m.*z_m + 3.*F(:,4).*z_m;

tmp31 = 0.5.*(x_m.^2 + y_m.^2 - 2.*z_m.^2).*F(:,5) - F(:,8).*y_m - F(:,10).*x_m + 0.5.*F(:,12) + 0.5.*F(:,13);
tmp32 = (x_m.^2 + y_m.^2 - 2.*z_m.^2).*F(:,4) - 2.*F(:,5).*y_m - 2.*F(:,6).*x_m + F(:,8) + F(:,9);
tmp33 = 0.5.*(x_m.^2 + y_m.^2 - 2.*z_m.^2).*F(:,3) - F(:,7).*x_m - F(:,9).*y_m + 0.5.*F(:,11) + 0.5.*F(:,12);
tmp34 = (x_m.^2 + y_m.^2 - 2.*z_m.^2).*F(:,2) - 2.*F(:,3).*x_m - 2.*F(:,6).*y_m + F(:,7) + F(:,10);
tmp35 = (x_m.^2 + y_m.^2 - 2.*z_m.^2).*F(:,1) - 2.*F(:,2).*x_m - 2.*F(:,4).*y_m + F(:,3) + F(:,5);

infl_new(1,1,:) = reshape(tmp11,1,1,[]);
infl_new(1,2,:) = reshape(tmp12,1,1,[]);
infl_new(1,3,:) = reshape(tmp13,1,1,[]);
infl_new(1,4,:) = reshape(tmp14,1,1,[]);
infl_new(1,5,:) = reshape(tmp15,1,1,[]);

infl_new(2,1,:) = reshape(tmp21,1,1,[]);
infl_new(2,2,:) = reshape(tmp22,1,1,[]);
infl_new(2,3,:) = reshape(tmp23,1,1,[]);
infl_new(2,4,:) = reshape(tmp24,1,1,[]);
infl_new(2,5,:) = reshape(tmp25,1,1,[]);

infl_new(3,1,:) = reshape(tmp31,1,1,[]);
infl_new(3,2,:) = reshape(tmp32,1,1,[]);
infl_new(3,3,:) = reshape(tmp33,1,1,[]);
infl_new(3,4,:) = reshape(tmp34,1,1,[]);
infl_new(3,5,:) = reshape(tmp35,1,1,[]);

%%
infl_loc = real(infl_new);
infl_loc(:,:,idx_flp) = -infl_loc(:,:,idx_flp);

idx_nan = find(reshape(sum(any(isnan(infl_new) | isinf(infl_new))),[],1) > 0);
if any(idx_nan)
   disp('Nan induction in fcnHDVEIND_DB'); 
end
% disp(['Inf or NaN induction: ', num2str(length(idx_nan))]);
infl_loc(:,:,idx_on_edge) = infl_loc(:,:,idx_on_edge).*0;
infl_loc(isnan(infl_loc) | isinf(infl_loc)) = 0;

% infl_loc(:,:,:) = infl_loc(:,:,:).*-1;
end