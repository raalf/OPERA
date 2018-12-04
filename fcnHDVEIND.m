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
C = (eta_3 - eta_2)./(xi_3 - xi_2);
D_LE = eta_2 - ((xi_2.*(eta_3 - eta_2))./(xi_3 - xi_2));
E = (eta_3 - eta_1)./(xi_3 - xi_1);
D_TE = eta_1 - ((xi_1.*(eta_3 - eta_1))./(xi_3 - xi_1));

le_eta = C.*x_m + D_LE;
te_eta = E.*x_m + D_TE;

margin_edge = 1e-10;
margin_on_element = 1e-10;
idx_on_element = y_m >= te_eta - margin_edge & y_m <= le_eta + margin_edge & x_m >= xi_1 - margin_edge & x_m <= xi_3 + margin_edge & abs(z_m) <= margin_on_element;
idx_on_edge =   (abs(y_m - te_eta) < margin_edge & (xi_1 - margin_edge <= x_m & x_m <= xi_3 + margin_edge)) | ...
                (abs(y_m - le_eta) < margin_edge & (xi_1 - margin_edge <= x_m & x_m <= xi_3 + margin_edge)) | ...
                (abs(x_m - xi_1) < margin_edge & (te_eta - margin_edge <= y_m & y_m <= le_eta + margin_edge)) | ...
                (abs(x_m - xi_3) < margin_edge & (te_eta - margin_edge <= y_m & y_m <= le_eta + margin_edge));
            
%% Calculating Influence
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

%% J_1
H_1_LE = fcnH_1(C, N_A, S_A, T_A, u_A, alpha, F1, F2);
H_1_TE = fcnH_1(E, N_B, S_B, T_B, u_B, alpha, F1, F2);
J_1 = -H_1_LE + H_1_TE;

%% J_2
q_a = y_m - 2.*C.*x_m - D_LE;
r_a = -x_m.*N_A;

q_b = y_m - 2.*E.*x_m - D_TE;
r_b = -x_m.*N_B;

H_2_LE = fcnH_2(S_A, T_A, u_A, F1, F2);
H_2_TE = fcnH_2(S_B, T_B, u_B, F1, F2);
J_2 = -(fcnH_1(q_a, r_a - C.*alpha, S_A, T_A, u_A, alpha, F1, F2) + C.*H_2_LE) + ...
      (fcnH_1(q_b, r_b - E.*alpha, S_B, T_B, u_B, alpha, F1, F2) + E.*H_2_TE);
J_2 = -J_2;
%% J_3
H_6_LE = fcnH_6(S_A, T_A, u_A, F1, F2);
H_6_TE = fcnH_6(S_B, T_B, u_B, F1, F2);

J_3 = -(fcnH_1((C.*x_m.^2 - C.*alpha + 2.*r_a), (-r_a.*x_m - (-2.*C.*x_m + N_A).*alpha), S_A, T_A, u_A, alpha, F1, F2) + (-2.*C.*x_m + N_A).*H_2_LE + C.*H_6_LE) + ...
      (fcnH_1((E.*x_m.^2 - E.*alpha + 2.*r_b), (-r_b.*x_m - (-2.*E.*x_m + N_B).*alpha), S_B, T_B, u_B, alpha, F1, F2) + (-2.*E.*x_m + N_B).*H_2_TE + E.*H_6_TE);  
  
%% J_4
J_4 = -(y_m.*H_1_LE + H_2_LE) + ...
      (y_m.*H_1_TE + H_2_TE);

%% J_5
J_5 = -((-y_m.^2).*H_1_LE + (N_A - 2.*y_m).*H_2_LE + C.*H_6_LE - fcnH_7(C, N_A, alpha, F1, F2)) + ...
       ((-y_m.^2).*H_1_TE + (N_B - 2.*y_m).*H_2_TE + E.*H_6_TE - fcnH_7(E, N_B, alpha, F1, F2));
J_5 = -J_5;
%% J_6
J_6 = -((C.*y_m - x_m).*H_2_LE + H_6_LE + fcnH_1((-C.*x_m.*y_m + N_A.*y_m), (-C.*alpha.*y_m - N_A.*x_m.*y_m), S_A, T_A, u_A, alpha, F1, F2)) + ...
       ((E.*y_m - x_m).*H_2_TE + H_6_TE + fcnH_1((-E.*x_m.*y_m + N_B.*y_m), (-E.*alpha.*y_m - N_B.*x_m.*y_m), S_B, T_B, u_B, alpha, F1, F2));
J_6 = -J_6;
%%
% J_1 = J_1.*0;
% J_2 = J_2.*0;
% J_3 = J_3.*0;
% J_4 = J_4.*0;
% J_5 = J_5.*0;
% J_6 = J_6.*0;
   
% J_1(idx_on_element) = J_1(idx_on_element).*0;
% J_2(idx_on_element) = J_1(idx_on_element).*0;
% J_3(idx_on_element) = J_1(idx_on_element).*0;
% J_4(idx_on_element) = J_1(idx_on_element).*0;
% J_5(idx_on_element) = J_1(idx_on_element).*0;
% J_6(idx_on_element) = J_1(idx_on_element).*0;

% % Whoopsie
% J_1(isnan(J_1) | isinf(J_1)) = 0;
% J_2(isnan(J_2) | isinf(J_2)) = 0;
% J_3(isnan(J_3) | isinf(J_3)) = 0;
% J_4(isnan(J_4) | isinf(J_4)) = 0;
% J_5(isnan(J_5) | isinf(J_5)) = 0;
% J_6(isnan(J_6) | isinf(J_6)) = 0;

% Compiling
infl_new = zeros(3,6,len);

infl_new(1,3,:) = reshape(-J_2.*z_m,1,1,[]);
infl_new(1,4,:) = reshape(J_1.*z_m,1,1,[]);
infl_new(1,5,:) = reshape(J_4.*z_m,1,1,[]);

infl_new(2,1,:) = reshape(J_4.*z_m,1,1,[]);
infl_new(2,2,:) = reshape(J_1.*z_m,1,1,[]);
infl_new(2,5,:) = reshape(J_2.*z_m,1,1,[]);

infl_new(3,1,:) = (reshape(-J_4.*y_m ,1,1,[]) + reshape(J_5,1,1,[]));
infl_new(3,2,:) = (reshape(-J_1.*y_m ,1,1,[]) + reshape(J_4,1,1,[]));
infl_new(3,3,:) = (reshape(J_2.*x_m,1,1,[]) + reshape(-J_3,1,1,[]));
infl_new(3,4,:) = (reshape(-J_1.*x_m,1,1,[]) + reshape(J_2,1,1,[]));
infl_new(3,5,:) = (reshape(-J_2.*y_m,1,1,[]) - reshape(J_4.*x_m,1,1,[]) +  reshape(2.*J_6,1,1,[]));

infl_loc = real(infl_new);

%% Transforming and Outputting
dvenum = reshape(repmat(dvenum,1,6,1)',[],1,1);
infl_tot = fcnSTARGLOB(reshape(permute(infl_loc,[2 3 1]),[],3,1), matROTANG(dvenum,:));

infl_tot(isnan(infl_tot)) = 0;
infl_tot(isinf(infl_tot)) = 0;

infl_glob = reshape(infl_tot',3,6,[]);

end