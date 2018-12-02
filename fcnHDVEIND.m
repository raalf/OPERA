function [infl_glob] = fcnHDVEIND(dvenum, dvetype, fpg, matPLEX, matROTANG, matCONTROL)
warning('on')
cutoff = 1e-14;

fpl = fcnGLOBSTAR(fpg - matCONTROL(dvenum,:), matROTANG(dvenum,:));

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
E = (eta_3 - eta_2)./(xi_3 - xi_2);
G = eta_2 - ((xi_2.*(eta_3 - eta_2))./(xi_3 - xi_2));
C = (eta_3 - eta_1)./(xi_3 - xi_1);
D = eta_1 - ((xi_1.*(eta_3 - eta_1))./(xi_3 - xi_1));

le_eta = E.*x_m + G;
te_eta = C.*x_m + D;

margin_edge = 1e-5;
margin_on_element = 1e-10;
idx_on_element = y_m >= te_eta - margin_edge & y_m <= le_eta + margin_edge & x_m >= xi_1 - margin_edge & x_m <= xi_3 + margin_edge & abs(z_m) <= margin_on_element;
idx_on_edge = abs(y_m - te_eta) < margin_edge | abs(y_m - le_eta) < margin_edge | abs(x_m - xi_1) < margin_edge | abs(x_m - xi_3) < margin_edge;

%% Calculating Influence
alpha = z_m.^2;
F1 = x_m;
F2 = x_m - xi_3;

% I_1
M_11 = C;
N_11 = y_m - C.*x_m - D;
S_11 = 1 + C.^2;
t_11 = 2.*C.*N_11;
u_11 = N_11.^2 + z_m.^2;
I_11 = fcnH_1(M_11, N_11, S_11, t_11, u_11, alpha, F1, F2);

M_12 = E;
N_12 = y_m - E.*x_m - G;
S_12 = 1 + E.^2;
t_12 = 2.*E.*N_12;
u_12 = N_12.^2 + z_m.^2;
I_12 = fcnH_1(M_12, N_12, S_12, t_12, u_12, alpha, F1, F2);

I_1 = I_11 + sign(F2 - F1).*I_12;

% I_2
P_21 = C;
q_21 = y_m - 2.*C.*x_m - D;
r_21 = -x_m.*(y_m - C.*x_m - D);
S_21 = 1 + C.^2;
t_21 = 2.*C.*(y_m - C.*x_m - D);
u_21 = (y_m - C.*x_m - D).^2 + z_m.^2;
I_21 = -P_21.*fcnH_2(S_21, t_21, u_21, F1, F2) - fcnH_1(q_21, (r_21 - P_21.*z_m.^2), S_21, t_21, u_21, alpha, F1, F2);

P_22 = E;
q_22 = y_m - 2.*E.*x_m - G;
r_22 = -x_m.*(y_m - E.*x_m - G);
S_22 = 1 + E.^2;
t_22 = 2.*E.*(y_m - E.*x_m - G);
u_22 = (y_m - E.*x_m - G).^2 + z_m.^2;
I_22 = P_22.*fcnH_2(S_22, t_22, u_22, F1, F2) + fcnH_1(q_22, (r_22 - P_22.*z_m.^2), S_22, t_22, u_22, alpha, F1, F2);

I_2 = I_21 + I_22;

% I_3
O_31 = C;
P_31 = y_m - 3.*C.*x_m - D;
q_31 = x_m.*(2.*D + 3.*C.*x_m - 2.*y_m);
r_31 = (x_m.^2).*(y_m - C.*x_m - D);
S_31 = 1 + C.^2;
t_31 = 2.*C.*(y_m - C.*x_m - D);
u_31 = (y_m - C.*x_m - D).^2 + z_m.^2;
M_31 = q_31 - O_31.*z_m.^2;
N_31 = r_31 - P_31.*z_m.^2;
I_31 = O_31.*fcnH_6(S_31, t_31, u_31, F1, F2) + P_31.*fcnH_2(S_31, t_31, u_31, F1, F2) + fcnH_1(M_31, N_31, S_31, t_31, u_31, alpha, F1, F2);

O_32 = E;
P_32 = y_m - 3.*E.*x_m - G;
q_32 = x_m.*(2.*G + 3.*E.*x_m - 2.*y_m);
r_32 = (x_m.^2).*(y_m - E.*x_m - G);
S_32 = 1 + E.^2;
t_32 = 2.*E.*(y_m - E.*x_m - G);
u_32 = (y_m - E.*x_m - G).^2 + z_m.^2;
M_32 = q_32 - O_32.*z_m.^2;
N_32 = r_32 - P_32.*z_m.^2;
I_32 = -(O_32.*fcnH_6(S_32, t_32, u_32, F1, F2) + P_32.*fcnH_2(S_32, t_32, u_32, F1, F2)) - fcnH_1(M_32, N_32, S_32, t_32, u_32, alpha, F1, F2);

I_3 = I_31 + I_32;

I_1(idx_on_element) = 0;
I_2(idx_on_element) = 0;
I_3(idx_on_element) = 0;

I_1(idx_on_edge) = 0;
I_2(idx_on_edge) = 0;
I_3(idx_on_edge) = 0;

% I_2(idx_on_element & ~idx_on_edge) = fcnH_8(x_m(idx_on_element & ~idx_on_edge), y_m(idx_on_element & ~idx_on_edge), C(idx_on_element & ~idx_on_edge), D(idx_on_element & ~idx_on_edge), E(idx_on_element & ~idx_on_edge), xi_1(idx_on_element & ~idx_on_edge), xi_3(idx_on_element & ~idx_on_edge));
% I_3(idx_on_element & ~idx_on_edge) = fcnH_9(x_m(idx_on_element & ~idx_on_edge), y_m(idx_on_element & ~idx_on_edge), C(idx_on_element & ~idx_on_edge), D(idx_on_element & ~idx_on_edge), E(idx_on_element & ~idx_on_edge), xi_1(idx_on_element & ~idx_on_edge), xi_3(idx_on_element & ~idx_on_edge));

%% Compiling
infl_new = zeros(3,3,len);

infl_new(2,1,:) = reshape(z_m.*I_1,1,1,[]);
infl_new(2,2,:) = reshape(z_m.*I_2,1,1,[]);

infl_new(3,1,:) = (reshape(-I_1.*y_m ,1,1,[]) + reshape(I_2,1,1,[]));
infl_new(3,2,:) = (reshape(-I_2.*y_m ,1,1,[]) + reshape(I_3,1,1,[]));

infl_loc = real(infl_new);

%% Transforming and Outputting
dvenum = reshape(repmat(dvenum,1,3,1)',[],1,1);
infl_tot = fcnSTARGLOB(reshape(permute(infl_loc,[2 3 1]),[],3,1), matROTANG(dvenum,:));

infl_tot(isnan(infl_tot)) = 0;
infl_tot(isinf(infl_tot)) = 0;

infl_glob = reshape(infl_tot',3,3,[]);

end