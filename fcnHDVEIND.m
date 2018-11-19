function [infl_glob] = fcnHDVEIND(dvenum, dvetype, fpg, matPLEX, matROTANG, matCENTER)
warning('on')
cutoff = 1e-14;

fpl = fcnGLOBSTAR(fpg - matCENTER(dvenum,:), matROTANG(dvenum,:));
len = size(fpl,1);

x_m = fpl(:,1);
y_m = fpl(:,2);
z_m = fpl(:,3);

%% Checking state of field point with relation to element surface
xi_1 = permute(matPLEX(1,1,dvenum),[3 2 1]);
xi_2 = permute(matPLEX(2,1,dvenum),[3 2 1]);
xi_3 = permute(matPLEX(3,1,dvenum),[3 2 1]);

eta_1 = permute(matPLEX(1,2,dvenum),[3 2 1]);
eta_2 = permute(matPLEX(2,2,dvenum),[3 2 1]);
eta_3 = permute(matPLEX(3,2,dvenum),[3 2 1]);

% Checking which elements are on the element
C = (eta_3 - eta_2)./(xi_3 - xi_2);
D_LE = eta_2 - ((xi_2.*(eta_3 - eta_2))./(xi_3 - xi_2));
E = (eta_3 - eta_1)./(xi_3 - xi_1);
D_TE = eta_1 - ((xi_1.*(eta_3 - eta_1))./(xi_3 - xi_1));

le_eta = C.*x_m + D_LE;
te_eta = E.*x_m + D_TE;

margin_edge = 1e-10;
margin_on_element = 1e-10;
idx_on_element = y_m >= te_eta - margin_edge & y_m <= le_eta + margin_edge & x_m >= xi_1 - margin_edge & x_m <= xi_3 + margin_edge & abs(z_m) <= margin_on_element;

%% Calculating Influence
tic

% J_1
alpha = z_m.^2;
N_A = -C.*x_m - D_LE + y_m;
S_A = C.^2 + 1;
T_A = 2.*C.*N_A;
u_A = N_A.^2 + alpha;

N_B = -E.*x_m - D_TE + y_m;
S_B = E.^2 + 1;
T_B = 2.*E.*N_B;
u_B = N_B.^2 + alpha;

F1 = (x_m - xi_1);
F2 = (x_m - xi_3);
J_1 = -fcnH_1(C, N_A, S_A, T_A, u_A, alpha, F1, F2) + fcnH_1(E, N_B, S_B, T_B, u_B, alpha, F1, F2);

q_a = y_m - 2.*C.*x_m - D_LE;
r_a = -x_m.*N_A;

q_b = y_m - 2.*E.*x_m - D_TE;
r_b = -x_m.*N_B;
% J_2 = -C.*fcnH_2(S_A, T_A, u_A, F1, F2) + E.*fcnH_2(S_B, T_B, u_B, F1, F2)
% J_2 = (-fcnH_1(q_a, r_a, S_A, T_A, u_A, alpha, F1, F2) - C.*fcnH_2(S_A, T_A, u_A, F1, F2)) + (fcnH_1(q_b, r_b, S_B, T_B, u_B, alpha, F1, F2) + E.*fcnH_2(S_B, T_B, u_B, F1, F2));
% J_1 = J_1.*0;






% Out of element plane
J_A0 = fcnJ_A0(x_m, y_m, z_m, xi_1, xi_3, C, D_LE);
J_A1 = fcnJ_A1(x_m, y_m, z_m, xi_1, xi_3, C, D_LE);
J_A2 = fcnJ_A2(x_m, y_m, z_m, xi_1, xi_3, C, D_LE);
J_A3 = fcnJ_A3(x_m, y_m, z_m, xi_1, xi_3, C, D_LE);

J_B0 = fcnJ_B0(x_m, y_m, z_m, xi_1, xi_3, E, D_TE);
J_B1 = fcnJ_B1(x_m, y_m, z_m, xi_1, xi_3, E, D_TE);
J_B2 = fcnJ_B2(x_m, y_m, z_m, xi_1, xi_3, E, D_TE);
J_B3 = fcnJ_B3(x_m, y_m, z_m, xi_1, xi_3, E, D_TE);

% J_1 = -C.*J_A1 - D_LE.*J_A0 + D_TE.*J_B0 + E.*J_B1 + J_A0.*y_m - J_B0.*y_m;
J_2 = -C.*J_A2 - D_LE.*J_A1 + D_TE.*J_B1 + E.*J_B2 + J_A1.*y_m - J_B1.*y_m;
J_3 = -C.*J_A3 - D_LE.*J_A2 + D_TE.*J_B2 + E.*J_B3 + J_A2.*y_m - J_B2.*y_m;
J_4 = -C.*J_A1.*y_m - D_LE.*J_A0.*y_m + D_TE.*J_B0.*y_m + E.*J_B1.*y_m + J_A0.*x_m.^2 + J_A0.*y_m.^2 + J_A0.*z_m.^2 - J_B0.*x_m.^2 - J_B0.*y_m.^2 - J_B0.*z_m.^2 - 2.*J_A1.*x_m + 2.*J_B1.*x_m + J_A2 - J_B2;
J_5 = (0.5.*((2.*J_A0-2.*J_B0).*y_m+2.*C.*J_A1+2.*D_LE.*J_A0-2.*D_TE.*J_B0-2.*E.*J_B1)).*x_m.^2+(0.5.*((-4.*J_A1+4.*J_B1).*y_m-4.*C.*J_A2+4.*E.*J_B2-4.*D_LE.*J_A1+4.*D_TE.*J_B1)).*x_m ...
    +(0.5.*(2.*J_A0-2.*J_B0)).*y_m.^3+(0.5.*(-2.*C.*J_A1-2.*D_LE.*J_A0+2.*D_TE.*J_B0+2.*E.*J_B1)).*y_m.^2+(0.5.*((2.*J_A0-2.*J_B0).*z_m.^2-2.*J_B2+2.*J_A2)).*y_m ...
    +(0.5.*(2.*C.*J_A1+2.*D_LE.*J_A0-2.*D_TE.*J_B0-2.*E.*J_B1)).*z_m.^2+J_A3.*C-J_B2.*D_TE-J_B3.*E+J_A2.*D_LE ...
    - (0.5.*fcnJ_D(x_m, y_m, z_m, xi_1, xi_3, C, D_LE)) + (0.5.*fcnJ_D(x_m, y_m, z_m, xi_1, xi_3, E, D_TE)) - (0.5.*fcnJ_C(x_m, y_m, z_m, xi_1, xi_3, C, D_LE)) + (0.5.*fcnJ_C(x_m, y_m, z_m, xi_1, xi_3, E, D_TE));
% J_5 = J_4.*0;
J_6 = -C.*J_A2.*y_m - D_LE.*J_A1.*y_m + D_TE.*J_B1.*y_m + E.*J_B2.*y_m + J_A1.*x_m.^2 + J_A1.*y_m.^2 + J_A1.*z_m.^2 - J_B1.*x_m.^2 - J_B1.*y_m.^2 - J_B1.*z_m.^2 - 2.*J_A2.*x_m + 2.*J_B2.*x_m + J_A3 - J_B3;


% % In plane of element (not on element)
% m_inf = 0;
% idx = fpl(:,3) == 0 & ~idx_on_element;
% idx_case1 = eta_3 <= eta_1;
% idx_case2 = eta_1 < eta_3 & eta_3 < eta_2;
% idx_case3 = eta_2 <= eta_3;
% J_1(idx) = fcnJ_1ip(E(idx), C(idx), D_LE(idx), D_TE(idx), x_m(idx), xi_1(idx), xi_3(idx), y_m(idx), cutoff, m_inf);
% J_2(idx) = fcnJ_2ip(E(idx), C(idx), D_LE(idx), D_TE(idx), x_m(idx), xi_1(idx), xi_3(idx), y_m(idx), cutoff, m_inf);
% J_3(idx) = fcnJ_3ip(E(idx), C(idx), D_LE(idx), D_TE(idx), x_m(idx), xi_1(idx), xi_3(idx), y_m(idx), cutoff, m_inf);
% J_4(idx) = fcnJ_4ip(E(idx), C(idx), D_LE(idx), D_TE(idx), x_m(idx), xi_1(idx), xi_3(idx), y_m(idx), cutoff, m_inf);
% J_5(idx & idx_case1) = fcnJ_51ip(x_m(idx & idx_case1), y_m(idx & idx_case1), xi_1(idx & idx_case1), eta_1(idx & idx_case1), eta_2(idx & idx_case1), eta_3(idx & idx_case1), C(idx & idx_case1), D_LE(idx & idx_case1), E(idx & idx_case1), D_TE(idx & idx_case1), cutoff, m_inf);
% J_5(idx & idx_case2) = fcnJ_52ip(x_m(idx & idx_case2), y_m(idx & idx_case2), xi_1(idx & idx_case2), eta_1(idx & idx_case2), eta_2(idx & idx_case2), eta_3(idx & idx_case2), C(idx & idx_case2), D_LE(idx & idx_case2), E(idx & idx_case2), D_TE(idx & idx_case2), cutoff, m_inf);
% J_5(idx & idx_case3) = fcnJ_53ip(x_m(idx & idx_case3), y_m(idx & idx_case3), xi_1(idx & idx_case3), eta_1(idx & idx_case3), eta_2(idx & idx_case3), eta_3(idx & idx_case3), C(idx & idx_case3), D_LE(idx & idx_case3), E(idx & idx_case3), D_TE(idx & idx_case3), cutoff, m_inf);
% J_6(idx) = fcnJ_6ip(E(idx), C(idx), D_LE(idx), D_TE(idx), x_m(idx), xi_1(idx), xi_3(idx), y_m(idx), cutoff, m_inf);

% % On element
% m_inf = 1;
% idx = fpl(:,3) == 0 & idx_on_element;
% idx_case1 = eta_3 <= eta_1;
% idx_case2 = eta_1 < eta_3 & eta_3 < eta_2;
% idx_case3 = eta_2 <= eta_3;
% J_1(idx) = fcnJ_1ip(E(idx), C(idx), D_LE(idx), D_TE(idx), x_m(idx), xi_1(idx), xi_3(idx), y_m(idx), cutoff, m_inf);
% J_2(idx) = fcnJ_2ip(E(idx), C(idx), D_LE(idx), D_TE(idx), x_m(idx), xi_1(idx), xi_3(idx), y_m(idx), cutoff, m_inf);
% J_3(idx) = fcnJ_3ip(E(idx), C(idx), D_LE(idx), D_TE(idx), x_m(idx), xi_1(idx), xi_3(idx), y_m(idx), cutoff, m_inf);
% J_4(idx) = fcnJ_4ip(E(idx), C(idx), D_LE(idx), D_TE(idx), x_m(idx), xi_1(idx), xi_3(idx), y_m(idx), cutoff, m_inf);
% J_5(idx & idx_case1) = fcnJ_51ip(x_m(idx & idx_case1), y_m(idx & idx_case1), xi_1(idx & idx_case1), eta_1(idx & idx_case1), eta_2(idx & idx_case1), eta_3(idx & idx_case1), C(idx & idx_case1), D_LE(idx & idx_case1), E(idx & idx_case1), D_TE(idx & idx_case1), cutoff, m_inf);
% J_5(idx & idx_case2) = fcnJ_52ip(x_m(idx & idx_case2), y_m(idx & idx_case2), xi_1(idx & idx_case2), eta_1(idx & idx_case2), eta_2(idx & idx_case2), eta_3(idx & idx_case2), C(idx & idx_case2), D_LE(idx & idx_case2), E(idx & idx_case2), D_TE(idx & idx_case2), cutoff, m_inf);
% J_5(idx & idx_case3) = fcnJ_53ip(x_m(idx & idx_case3), y_m(idx & idx_case3), xi_1(idx & idx_case3), eta_1(idx & idx_case3), eta_2(idx & idx_case3), eta_3(idx & idx_case3), C(idx & idx_case3), D_LE(idx & idx_case3), E(idx & idx_case3), D_TE(idx & idx_case3), cutoff, m_inf);
% J_6(idx) = fcnJ_6ip(E(idx), C(idx), D_LE(idx), D_TE(idx), x_m(idx), xi_1(idx), xi_3(idx), y_m(idx), cutoff, m_inf);

% Whoopsie
J_1(isnan(J_1) | isinf(J_1)) = 0;
J_2(isnan(J_2) | isinf(J_2)) = 0;
J_3(isnan(J_3) | isinf(J_3)) = 0;
J_4(isnan(J_4) | isinf(J_4)) = 0;
J_5(isnan(J_5) | isinf(J_5)) = 0;
J_6(isnan(J_6) | isinf(J_6)) = 0;

% Compiling
infl_new = zeros(3,6,len);

infl_new(1,3,:) = reshape(J_2.*z_m,1,1,[]);
infl_new(1,4,:) = reshape(J_1.*z_m,1,1,[]);
infl_new(1,5,:) = reshape(J_4.*z_m,1,1,[]);

infl_new(2,1,:) = reshape(J_4.*z_m,1,1,[]);
infl_new(2,2,:) = reshape(J_1.*z_m,1,1,[]);
infl_new(2,5,:) = reshape(J_2.*z_m,1,1,[]);

infl_new(3,1,:) = (reshape(-J_4.*y_m ,1,1,[]) + reshape(J_5,1,1,[]));
infl_new(3,2,:) = (reshape(-J_1.*y_m ,1,1,[]) + reshape(J_4,1,1,[]));
infl_new(3,3,:) = (reshape(-J_2.*x_m,1,1,[]) + reshape(J_3,1,1,[]));
infl_new(3,4,:) = (reshape(-J_1.*x_m,1,1,[]) + reshape(J_2,1,1,[]));
infl_new(3,5,:) = (reshape(-J_2.*y_m,1,1,[]) - reshape(J_4.*x_m,1,1,[]) +  reshape(2.*J_6,1,1,[]));

infl_loc = real(infl_new);
toc

%% Transforming and Outputting
dvenum = reshape(repmat(dvenum,1,6,1)',[],1,1);
infl_tot = fcnSTARGLOB(reshape(permute(infl_loc,[2 3 1]),[],3,1), matROTANG(dvenum,:));

infl_tot(isnan(infl_tot)) = 0;
infl_tot(isinf(infl_tot)) = 0;

infl_glob = reshape(infl_tot',3,6,[]);

end

function J_1ip = fcnJ_1ip(E, C, D_LE, D_TE, x_m, xi_1, xi_3, y_m, cutoff, m_inf)
J_1ip = E.*nan;
end

function J_2ip = fcnJ_2ip(E, C, D_LE, D_TE, x_m, xi_1, xi_3, y_m, cutoff, m_inf)
J_2ip = E.*nan;
end

function J_3ip = fcnJ_3ip(E, C, D_LE, D_TE, x_m, xi_1, xi_3, y_m, cutoff, m_inf)
J_3ip = E.*nan;
end

function J_4ip = fcnJ_4ip(E, C, D_LE, D_TE, x_m, xi_1, xi_3, y_m, cutoff, m_inf)
J_4ip = E.*nan;
end

function J_5ip = fcnJ_5ip(E, C, D_LE, D_TE, x_m, xi_1, xi_3, y_m, cutoff, m_inf)
J_5ip = E.*nan;
end

function J_6ip = fcnJ_6ip(E, C, D_LE, D_TE, x_m, xi_1, xi_3, y_m, cutoff, m_inf)
J_6ip = E.*nan;
end

function out = fcnPIECEWISE(cond,val1,val2)
out = val2;
out(cond) = val1(cond);
end

function out = fcnAND(cond1, cond2)
out = cond1 & cond2;
end
