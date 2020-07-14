function [infl_loc] = fcnBOUNDIND(dvenum, dvetype, fpg, matPLEX, matROTANG, matCONTROL, vecDVESDFLIP, vecBI, ztol, GPU)
warning('on')
% tol = 1e-10;

if nargin < 9 || GPU == false
    GPU = false;
end

fpl = fcnGLOBSTAR(fpg - matCONTROL(dvenum,:), matROTANG(dvenum,:), GPU);

x_m = fpl(:,1);
y_m = fpl(:,2);
z_m = fpl(:,3);

idx_afx = abs(z_m) < ztol;
z_m_orig = z_m(idx_afx);
sgn = sign(z_m_orig);
sgn(sgn == 0) = 1;
z_m(idx_afx) = sgn.*ztol;

%% Checking state of field point with relation to element surface
margin_edge = 1e-10;

xi_1 = permute(matPLEX(1,1,dvenum),[3 2 1]);
xi_2 = permute(matPLEX(2,1,dvenum),[3 2 1]);
xi_3 = permute(matPLEX(3,1,dvenum),[3 2 1]);

eta_1 = permute(matPLEX(1,2,dvenum),[3 2 1]);
eta_2 = permute(matPLEX(2,2,dvenum),[3 2 1]);
eta_3 = permute(matPLEX(3,2,dvenum),[3 2 1]);

idx_flp = xi_3 < xi_1; % Flipping influence of elements that need a good flippin
% idx_flp(vecDVESDFLIP) = ~idx_flp(vecDVESDFLIP);
if any(abs(xi_2 - xi_3) < margin_edge & abs(xi_1 - xi_2) > margin_edge)
    disp('Issue in element orientation in HDVEIND.');
end

% Checking which elements are on the element
C = (eta_3 - eta_2)./(xi_3 - xi_2);
D_LE = eta_2 - ((xi_2.*(eta_3 - eta_2))./(xi_3 - xi_2));
E = (eta_3 - eta_1)./(xi_3 - xi_1);
D_TE = eta_1 - ((xi_1.*(eta_3 - eta_1))./(xi_3 - xi_1));

%% Bound induction

if any(vecBI)
    t1 = (eta_2(vecBI) .^ 2);
    t2 = (y_m(vecBI) .* eta_2(vecBI));
    t4 = (x_m(vecBI) .^ 2);
    t7 = (xi_2(vecBI) .^ 2);
    t8 = (y_m(vecBI) .^ 2);
    t9 = (z_m(vecBI) .^ 2);
    t11 = sqrt((-2 .* x_m(vecBI) .* xi_2(vecBI) + t1 - 2 .* t2 + t4 + t7 + t8 + t9));
    t14 = (eta_3(vecBI) .^ 2);
    t15 = (y_m(vecBI) .* eta_3(vecBI));
    t17 = x_m(vecBI) .* xi_3(vecBI);
    t19 = xi_3(vecBI) .^ 2;
    t21 = sqrt((t14 - 2 .* t15 + t4 - 2 .* t17 + t19 + t8 + t9));
    t23 = (C(vecBI) .^ 2);
    t34 = (D_LE(vecBI) .^ 2);
    t41 = (eta_2(vecBI) - eta_3(vecBI));
    t42 = (y_m(vecBI) - eta_3(vecBI));
    t49 = ((y_m(vecBI) - eta_2(vecBI)) .* t41);
    t56 = 12 .* t15;
    t57 = 6 .* t8;
    t58 = 6 .* t9;
    t78 = (t41 .^ 2);
    tmp11 = 0.1e1 ./ (t7 .* (6 .* t14 - t56 + t57 + t58) + xi_2(vecBI) .* (xi_3(vecBI) .* (12 .* eta_2(vecBI) .* t42 + t56 - 12 .* t8 - 12 .* t9) - 12 .* t42 .* t41 .* x_m(vecBI)) + t19 .* (6 .* t1 - 12 .* t2 + t57 + t58) + 12 .* t17 .* t49 + 6 .* (t4 + t9) .* t78) .* t41 .* (-xi_3(vecBI) + xi_2(vecBI)) .* (t11 .* ((x_m(vecBI) - xi_3(vecBI)) .* xi_2(vecBI) + t19 - t17 + t42 .* t41) - (-t7 + (xi_3(vecBI) + x_m(vecBI)) .* xi_2(vecBI) - t17 + t49) .* t21) .* (t7 .* t23 + C(vecBI) .* (C(vecBI) .* xi_3(vecBI) + 3 .* D_LE(vecBI)) .* xi_2(vecBI) + t19 .* t23 + 3 .* C(vecBI) .* D_LE(vecBI) .* xi_3(vecBI) + 3 .* t34) ./ t21 .* z_m(vecBI) ./ t11;
    t1 = (eta_2(vecBI) .^ 2);
    t2 = (y_m(vecBI) .* eta_2(vecBI));
    t4 = (x_m(vecBI) .^ 2);
    t7 = (xi_2(vecBI) .^ 2);
    t8 = (y_m(vecBI) .^ 2);
    t9 = (z_m(vecBI) .^ 2);
    t11 = sqrt((-2 .* x_m(vecBI) .* xi_2(vecBI) + t1 - 2 .* t2 + t4 + t7 + t8 + t9));
    t14 = (eta_3(vecBI) .^ 2);
    t15 = (y_m(vecBI) .* eta_3(vecBI));
    t17 = x_m(vecBI) .* xi_3(vecBI);
    t19 = xi_3(vecBI) .^ 2;
    t21 = sqrt((t14 - 2 .* t15 + t4 - 2 .* t17 + t19 + t8 + t9));
    t25 = (eta_2(vecBI) - eta_3(vecBI));
    t26 = (y_m(vecBI) - eta_3(vecBI));
    t33 = ((y_m(vecBI) - eta_2(vecBI)) .* t25);
    t46 = 4 .* t15;
    t47 = 2 .* t8;
    t48 = 2 .* t9;
    t68 = (t25 .^ 2);
    tmp12 = 0.1e1 ./ (t7 .* (2 .* t14 - t46 + t47 + t48) + xi_2(vecBI) .* (xi_3(vecBI) .* (4 .* eta_2(vecBI) .* t26 + t46 - 4 .* t8 - 4 .* t9) - 4 .* t26 .* t25 .* x_m(vecBI)) + t19 .* (2 .* t1 - 4 .* t2 + t47 + t48) + 4 .* t17 .* t33 + 2 .* (t4 + t9) .* t68) .* t25 .* (C(vecBI) .* xi_2(vecBI) + C(vecBI) .* xi_3(vecBI) + 2 .* D_LE(vecBI)) .* (-xi_3(vecBI) + xi_2(vecBI)) .* (t11 .* ((x_m(vecBI) - xi_3(vecBI)) .* xi_2(vecBI) + t19 - t17 + t26 .* t25) - (-t7 + (xi_3(vecBI) + x_m(vecBI)) .* xi_2(vecBI) - t17 + t33) .* t21) ./ t21 .* z_m(vecBI) ./ t11;
    t1 = (-xi_3(vecBI) + xi_2(vecBI));
    t2 = (eta_3(vecBI) .^ 2);
    t3 = (y_m(vecBI) .* eta_3(vecBI));
    t5 = (x_m(vecBI) .^ 2);
    t6 = (x_m(vecBI) .* xi_3(vecBI));
    t8 = (xi_3(vecBI) .^ 2);
    t9 = (y_m(vecBI) .^ 2);
    t10 = (z_m(vecBI) .^ 2);
    t12 = sqrt((t2 - 2 .* t3 + t5 - 2 .* t6 + t8 + t9 + t10));
    t16 = (eta_2(vecBI) .^ 2);
    t19 = (x_m(vecBI) .* xi_2(vecBI));
    t21 = (xi_2(vecBI) .^ 2);
    t23 = sqrt((-2 .* y_m(vecBI) .* eta_2(vecBI) + t10 + t16 - 2 .* t19 + t21 + t5 + t9));
    t29 = (x_m(vecBI) - xi_3(vecBI));
    t35 = (x_m(vecBI) - xi_2(vecBI));
    t43 = 6 .* t5;
    t44 = 12 .* t6;
    t46 = 6 .* t10;
    t55 = t1 .* y_m(vecBI);
    t67 = t1 .^ 2;
    tmp13 = -0.1e1 ./ (t16 .* (t43 - t44 + 6 .* t8 + t46) + eta_2(vecBI) .* (eta_3(vecBI) .* (12 .* xi_2(vecBI) .* t29 - 12 .* t10 + t44 - 12 .* t5) - 12 .* t29 .* t55) + t2 .* (t43 - 12 .* t19 + 6 .* t21 + t46) + 12 .* eta_3(vecBI) .* t35 .* t55 + 6 .* (t9 + t10) .* t67) .* (xi_3(vecBI) .* xi_2(vecBI) + t21 + t8) .* (t23 .* ((-y_m(vecBI) + eta_3(vecBI)) .* eta_2(vecBI) - t2 + t3 - t29 .* t1) + t12 .* (-t16 + (y_m(vecBI) + eta_3(vecBI)) .* eta_2(vecBI) - t3 + t35 .* t1)) .* z_m(vecBI) ./ t23 .* (eta_2(vecBI) - eta_3(vecBI)) ./ t12 .* t1;
    t1 = (-xi_3(vecBI) + xi_2(vecBI));
    t2 = (eta_3(vecBI) .^ 2);
    t3 = (y_m(vecBI) .* eta_3(vecBI));
    t5 = (x_m(vecBI) .^ 2);
    t6 = (x_m(vecBI) .* xi_3(vecBI));
    t8 = (xi_3(vecBI) .^ 2);
    t9 = (y_m(vecBI) .^ 2);
    t10 = (z_m(vecBI) .^ 2);
    t12 = sqrt((t2 - 2 .* t3 + t5 - 2 .* t6 + t8 + t9 + t10));
    t19 = (eta_2(vecBI) .^ 2);
    t22 = (x_m(vecBI) .* xi_2(vecBI));
    t24 = (xi_2(vecBI) .^ 2);
    t26 = sqrt((-2 .* y_m(vecBI) .* eta_2(vecBI) + t10 + t19 - 2 .* t22 + t24 + t5 + t9));
    t31 = (x_m(vecBI) - xi_3(vecBI));
    t37 = (x_m(vecBI) - xi_2(vecBI));
    t42 = 2 .* t5;
    t43 = 4 .* t6;
    t45 = 2 .* t10;
    t54 = t1 .* y_m(vecBI);
    t66 = t1 .^ 2;
    tmp14 = -0.1e1 ./ (t19 .* (t42 - t43 + 2 .* t8 + t45) + eta_2(vecBI) .* (eta_3(vecBI) .* (4 .* xi_2(vecBI) .* t31 - 4 .* t10 + t43 - 4 .* t5) - 4 .* t31 .* t54) + t2 .* (t42 - 4 .* t22 + 2 .* t24 + t45) + 4 .* eta_3(vecBI) .* t37 .* t54 + 2 .* (t9 + t10) .* t66) .* (t26 .* ((-y_m(vecBI) + eta_3(vecBI)) .* eta_2(vecBI) - t2 + t3 - t31 .* t1) + t12 .* (-t19 + (y_m(vecBI) + eta_3(vecBI)) .* eta_2(vecBI) - t3 + t37 .* t1)) .* z_m(vecBI) ./ t26 .* (xi_2(vecBI) + xi_3(vecBI)) .* (eta_2(vecBI) - eta_3(vecBI)) ./ t12 .* t1;
    t1 = (eta_2(vecBI) .^ 2);
    t2 = (y_m(vecBI) .* eta_2(vecBI));
    t4 = (x_m(vecBI) .^ 2);
    t7 = (xi_2(vecBI) .^ 2);
    t8 = (y_m(vecBI) .^ 2);
    t9 = (z_m(vecBI) .^ 2);
    t11 = sqrt((-2 .* x_m(vecBI) .* xi_2(vecBI) + t1 - 2 .* t2 + t4 + t7 + t8 + t9));
    t14 = (eta_3(vecBI) .^ 2);
    t15 = (y_m(vecBI) .* eta_3(vecBI));
    t17 = x_m(vecBI) .* xi_3(vecBI);
    t19 = xi_3(vecBI) .^ 2;
    t21 = sqrt((t14 - 2 .* t15 + t4 - 2 .* t17 + t19 + t8 + t9));
    t25 = (eta_2(vecBI) - eta_3(vecBI));
    t26 = (y_m(vecBI) - eta_3(vecBI));
    t33 = ((y_m(vecBI) - eta_2(vecBI)) .* t25);
    t43 = (C(vecBI) .* xi_3(vecBI)) + 0.3e1 ./ 0.2e1 .* D_LE(vecBI);
    t49 = 12 .* t15;
    t50 = 6 .* t8;
    t51 = 6 .* t9;
    t71 = (t25 .^ 2);
    tmp15 = 0.2e1 ./ (t7 .* (6 .* t14 - t49 + t50 + t51) + xi_2(vecBI) .* (xi_3(vecBI) .* (12 .* eta_2(vecBI) .* t26 + t49 - 12 .* t8 - 12 .* t9) - 12 .* t26 .* t25 .* x_m(vecBI)) + t19 .* (6 .* t1 - 12 .* t2 + t50 + t51) + 12 .* t17 .* t33 + 6 .* (t4 + t9) .* t71) .* t25 .* ((t7 .* C(vecBI)) + xi_2(vecBI) .* t43 + t43 .* xi_3(vecBI)) .* (-xi_3(vecBI) + xi_2(vecBI)) .* (t11 .* ((x_m(vecBI) - xi_3(vecBI)) .* xi_2(vecBI) + t19 - t17 + t26 .* t25) - (-t7 + (xi_3(vecBI) + x_m(vecBI)) .* xi_2(vecBI) - t17 + t33) .* t21) ./ t21 .* z_m(vecBI) ./ t11;
    t1 = (-xi_3(vecBI) + xi_2(vecBI));
    t2 = (eta_3(vecBI) .^ 2);
    t3 = (y_m(vecBI) .* eta_3(vecBI));
    t5 = (x_m(vecBI) .^ 2);
    t7 = 2 .* x_m(vecBI) .* xi_3(vecBI);
    t8 = xi_3(vecBI) .^ 2;
    t9 = (y_m(vecBI) .^ 2);
    t10 = (z_m(vecBI) .^ 2);
    t12 = sqrt((t2 - 2 .* t3 + t5 - t7 + t8 + t9 + t10));
    t17 = (eta_2(vecBI) .^ 2);
    t21 = 2 .* x_m(vecBI) .* xi_2(vecBI);
    t22 = xi_2(vecBI) .^ 2;
    t24 = sqrt((-2 .* y_m(vecBI) .* eta_2(vecBI) + t10 + t17 - t21 + t22 + t5 + t9));
    t29 = x_m(vecBI) - xi_3(vecBI);
    t35 = x_m(vecBI) - xi_2(vecBI);
    t48 = t1 .* y_m(vecBI);
    t58 = t1 .^ 2;
    tmp16 = -0.1e1 ./ (t17 .* (t5 - t7 + t8 + t10) + eta_2(vecBI) .* (eta_3(vecBI) .* (2 .* xi_2(vecBI) .* t29 - 2 .* t10 - 2 .* t5 + t7) - 2 .* t29 .* t48) + t2 .* (t5 - t21 + t22 + t10) + 2 .* eta_3(vecBI) .* t35 .* t48 + (t9 + t10) .* t58) .* (t24 .* ((-y_m(vecBI) + eta_3(vecBI)) .* eta_2(vecBI) - t2 + t3 - t29 .* t1) + t12 .* (-t17 + (y_m(vecBI) + eta_3(vecBI)) .* eta_2(vecBI) - t3 + t35 .* t1)) .* z_m(vecBI) ./ t24 .* (eta_2(vecBI) - eta_3(vecBI)) ./ t12 .* t1;
    t1 = (eta_2(vecBI) .^ 2);
    t2 = (y_m(vecBI) .* eta_2(vecBI));
    t4 = (x_m(vecBI) .^ 2);
    t7 = (xi_2(vecBI) .^ 2);
    t8 = (y_m(vecBI) .^ 2);
    t9 = (z_m(vecBI) .^ 2);
    t11 = sqrt((-2 .* x_m(vecBI) .* xi_2(vecBI) + t1 - 2 .* t2 + t4 + t7 + t8 + t9));
    t14 = (eta_3(vecBI) .^ 2);
    t15 = (y_m(vecBI) .* eta_3(vecBI));
    t17 = x_m(vecBI) .* xi_3(vecBI);
    t19 = xi_3(vecBI) .^ 2;
    t21 = sqrt((t14 - 2 .* t15 + t4 - 2 .* t17 + t19 + t8 + t9));
    t24 = (C(vecBI) .^ 2);
    t35 = (D_LE(vecBI) .^ 2);
    t40 = (eta_2(vecBI) - eta_3(vecBI));
    t41 = (y_m(vecBI) - eta_3(vecBI));
    t48 = ((y_m(vecBI) - eta_2(vecBI)) .* t40);
    t54 = (-xi_3(vecBI) + xi_2(vecBI)) .^ 2;
    t56 = 12 .* t15;
    t57 = 6 .* t8;
    t58 = 6 .* t9;
    t78 = (t40 .^ 2);
    tmp21 = -0.1e1 ./ (t7 .* (6 .* t14 - t56 + t57 + t58) + xi_2(vecBI) .* (xi_3(vecBI) .* (12 .* eta_2(vecBI) .* t41 + t56 - 12 .* t8 - 12 .* t9) - 12 .* t41 .* t40 .* x_m(vecBI)) + t19 .* (6 .* t1 - 12 .* t2 + t57 + t58) + 12 .* t17 .* t48 + 6 .* (t4 + t9) .* t78) .* t54 .* (t11 .* ((x_m(vecBI) - xi_3(vecBI)) .* xi_2(vecBI) + t19 - t17 + t41 .* t40) - (-t7 + (xi_3(vecBI) + x_m(vecBI)) .* xi_2(vecBI) - t17 + t48) .* t21) .* (t7 .* t24 + C(vecBI) .* (C(vecBI) .* xi_3(vecBI) + 3 .* D_LE(vecBI)) .* xi_2(vecBI) + t19 .* t24 + 3 .* C(vecBI) .* D_LE(vecBI) .* xi_3(vecBI) + 3 .* t35) ./ t21 .* z_m(vecBI) ./ t11;
    t1 = (eta_2(vecBI) .^ 2);
    t2 = (y_m(vecBI) .* eta_2(vecBI));
    t4 = (x_m(vecBI) .^ 2);
    t7 = (xi_2(vecBI) .^ 2);
    t8 = (y_m(vecBI) .^ 2);
    t9 = (z_m(vecBI) .^ 2);
    t11 = sqrt((-2 .* x_m(vecBI) .* xi_2(vecBI) + t1 - 2 .* t2 + t4 + t7 + t8 + t9));
    t14 = (eta_3(vecBI) .^ 2);
    t15 = (y_m(vecBI) .* eta_3(vecBI));
    t17 = x_m(vecBI) .* xi_3(vecBI);
    t19 = xi_3(vecBI) .^ 2;
    t21 = sqrt((t14 - 2 .* t15 + t4 - 2 .* t17 + t19 + t8 + t9));
    t26 = (eta_2(vecBI) - eta_3(vecBI));
    t27 = (y_m(vecBI) - eta_3(vecBI));
    t34 = ((y_m(vecBI) - eta_2(vecBI)) .* t26);
    t44 = (-xi_3(vecBI) + xi_2(vecBI)) .^ 2;
    t46 = 4 .* t15;
    t47 = 2 .* t8;
    t48 = 2 .* t9;
    t68 = (t26 .^ 2);
    tmp22 = -0.1e1 ./ (t7 .* (2 .* t14 - t46 + t47 + t48) + xi_2(vecBI) .* (xi_3(vecBI) .* (4 .* eta_2(vecBI) .* t27 + t46 - 4 .* t8 - 4 .* t9) - 4 .* t27 .* t26 .* x_m(vecBI)) + t19 .* (2 .* t1 - 4 .* t2 + t47 + t48) + 4 .* t17 .* t34 + 2 .* (t4 + t9) .* t68) .* t44 .* (C(vecBI) .* xi_2(vecBI) + C(vecBI) .* xi_3(vecBI) + 2 .* D_LE(vecBI)) .* (t11 .* ((x_m(vecBI) - xi_3(vecBI)) .* xi_2(vecBI) + t19 - t17 + t27 .* t26) - (-t7 + (xi_3(vecBI) + x_m(vecBI)) .* xi_2(vecBI) - t17 + t34) .* t21) ./ t21 .* z_m(vecBI) ./ t11;
    t2 = (-xi_3(vecBI) + xi_2(vecBI)) .^ 2;
    t3 = (eta_3(vecBI) .^ 2);
    t4 = (y_m(vecBI) .* eta_3(vecBI));
    t6 = (x_m(vecBI) .^ 2);
    t7 = (x_m(vecBI) .* xi_3(vecBI));
    t9 = (xi_3(vecBI) .^ 2);
    t10 = (y_m(vecBI) .^ 2);
    t11 = (z_m(vecBI) .^ 2);
    t13 = sqrt((t3 - 2 .* t4 + t6 - 2 .* t7 + t9 + t10 + t11));
    t16 = (eta_2(vecBI) .^ 2);
    t17 = (y_m(vecBI) .* eta_2(vecBI));
    t21 = (xi_2(vecBI) .^ 2);
    t23 = sqrt((-2 .* x_m(vecBI) .* xi_2(vecBI) + t10 + t11 + t16 - 2 .* t17 + t21 + t6));
    t28 = (eta_2(vecBI) - eta_3(vecBI));
    t29 = (y_m(vecBI) - eta_3(vecBI));
    t36 = ((y_m(vecBI) - eta_2(vecBI)) .* t28);
    t44 = 12 .* t4;
    t45 = 6 .* t10;
    t46 = 6 .* t11;
    t66 = (t28 .^ 2);
    tmp23 = 0.1e1 ./ (t21 .* (6 .* t3 - t44 + t45 + t46) + xi_2(vecBI) .* (xi_3(vecBI) .* (12 .* eta_2(vecBI) .* t29 - 12 .* t10 - 12 .* t11 + t44) - 12 .* t29 .* t28 .* x_m(vecBI)) + t9 .* (6 .* t16 - 12 .* t17 + t45 + t46) + 12 .* t7 .* t36 + 6 .* (t6 + t11) .* t66) .* (xi_2(vecBI) .* xi_3(vecBI) + t21 + t9) .* (t23 .* ((xi_3(vecBI) - x_m(vecBI)) .* xi_2(vecBI) - t9 + t7 - t29 .* t28) + (-t21 + (xi_3(vecBI) + x_m(vecBI)) .* xi_2(vecBI) - t7 + t36) .* t13) .* z_m(vecBI) ./ t23 ./ t13 .* t2;
    t2 = (-xi_3(vecBI) + xi_2(vecBI)) .^ 2;
    t3 = (eta_3(vecBI) .^ 2);
    t4 = (y_m(vecBI) .* eta_3(vecBI));
    t6 = (x_m(vecBI) .^ 2);
    t7 = (x_m(vecBI) .* xi_3(vecBI));
    t9 = (xi_3(vecBI) .^ 2);
    t10 = (y_m(vecBI) .^ 2);
    t11 = (z_m(vecBI) .^ 2);
    t13 = sqrt((t3 - 2 .* t4 + t6 - 2 .* t7 + t9 + t10 + t11));
    t18 = (eta_2(vecBI) .^ 2);
    t19 = (y_m(vecBI) .* eta_2(vecBI));
    t23 = (xi_2(vecBI) .^ 2);
    t25 = sqrt((-2 .* x_m(vecBI) .* xi_2(vecBI) + t10 + t11 + t18 - 2 .* t19 + t23 + t6));
    t30 = (eta_2(vecBI) - eta_3(vecBI));
    t31 = (y_m(vecBI) - eta_3(vecBI));
    t38 = ((y_m(vecBI) - eta_2(vecBI)) .* t30);
    t43 = 4 .* t4;
    t44 = 2 .* t10;
    t45 = 2 .* t11;
    t65 = (t30 .^ 2);
    tmp24 = 0.1e1 ./ (t23 .* (2 .* t3 - t43 + t44 + t45) + xi_2(vecBI) .* (xi_3(vecBI) .* (4 .* eta_2(vecBI) .* t31 - 4 .* t10 - 4 .* t11 + t43) - 4 .* t31 .* t30 .* x_m(vecBI)) + t9 .* (2 .* t18 - 4 .* t19 + t44 + t45) + 4 .* t7 .* t38 + 2 .* (t6 + t11) .* t65) .* (t25 .* ((xi_3(vecBI) - x_m(vecBI)) .* xi_2(vecBI) - t9 + t7 - t31 .* t30) + (-t23 + (xi_3(vecBI) + x_m(vecBI)) .* xi_2(vecBI) - t7 + t38) .* t13) .* z_m(vecBI) ./ t25 .* (xi_2(vecBI) + xi_3(vecBI)) ./ t13 .* t2;
    t1 = (eta_2(vecBI) .^ 2);
    t2 = (y_m(vecBI) .* eta_2(vecBI));
    t4 = (x_m(vecBI) .^ 2);
    t7 = (xi_2(vecBI) .^ 2);
    t8 = (y_m(vecBI) .^ 2);
    t9 = (z_m(vecBI) .^ 2);
    t11 = sqrt((-2 .* x_m(vecBI) .* xi_2(vecBI) + t1 - 2 .* t2 + t4 + t7 + t8 + t9));
    t14 = (eta_3(vecBI) .^ 2);
    t15 = (y_m(vecBI) .* eta_3(vecBI));
    t17 = x_m(vecBI) .* xi_3(vecBI);
    t19 = xi_3(vecBI) .^ 2;
    t21 = sqrt((t14 - 2 .* t15 + t4 - 2 .* t17 + t19 + t8 + t9));
    t26 = (eta_2(vecBI) - eta_3(vecBI));
    t27 = (y_m(vecBI) - eta_3(vecBI));
    t34 = ((y_m(vecBI) - eta_2(vecBI)) .* t26);
    t39 = (-xi_3(vecBI) + xi_2(vecBI)) .^ 2;
    t44 = (C(vecBI) .* xi_3(vecBI)) + 0.3e1 ./ 0.2e1 .* D_LE(vecBI);
    t49 = 12 .* t15;
    t50 = 6 .* t8;
    t51 = 6 .* t9;
    t71 = (t26 .^ 2);
    tmp25 = -0.2e1 ./ (t7 .* (6 .* t14 - t49 + t50 + t51) + xi_2(vecBI) .* (xi_3(vecBI) .* (12 .* eta_2(vecBI) .* t27 + t49 - 12 .* t8 - 12 .* t9) - 12 .* t27 .* t26 .* x_m(vecBI)) + t19 .* (6 .* t1 - 12 .* t2 + t50 + t51) + 12 .* t17 .* t34 + 6 .* (t4 + t9) .* t71) .* ((t7 .* C(vecBI)) + xi_2(vecBI) .* t44 + t44 .* xi_3(vecBI)) .* t39 .* (t11 .* ((x_m(vecBI) - xi_3(vecBI)) .* xi_2(vecBI) + t19 - t17 + t27 .* t26) - (-t7 + (xi_3(vecBI) + x_m(vecBI)) .* xi_2(vecBI) - t17 + t34) .* t21) ./ t21 .* z_m(vecBI) ./ t11;
    t2 = (-xi_3(vecBI) + xi_2(vecBI)) .^ 2;
    t3 = (eta_3(vecBI) .^ 2);
    t5 = 2 .* y_m(vecBI) .* eta_3(vecBI);
    t6 = (x_m(vecBI) .^ 2);
    t7 = (x_m(vecBI) .* xi_3(vecBI));
    t9 = (xi_3(vecBI) .^ 2);
    t10 = y_m(vecBI) .^ 2;
    t11 = (z_m(vecBI) .^ 2);
    t13 = sqrt((t3 - t5 + t6 - 2 .* t7 + t9 + t10 + t11));
    t16 = (eta_2(vecBI) .^ 2);
    t18 = 2 .* y_m(vecBI) .* eta_2(vecBI);
    t21 = (xi_2(vecBI) .^ 2);
    t23 = sqrt((-2 .* x_m(vecBI) .* xi_2(vecBI) + t10 + t11 + t16 - t18 + t21 + t6));
    t28 = eta_2(vecBI) - eta_3(vecBI);
    t29 = y_m(vecBI) - eta_3(vecBI);
    t36 = (y_m(vecBI) - eta_2(vecBI)) .* t28;
    t58 = t28 .^ 2;
    tmp26 = 0.1e1 ./ (t21 .* (t3 - t5 + t10 + t11) + xi_2(vecBI) .* (xi_3(vecBI) .* (2 .* eta_2(vecBI) .* t29 - 2 .* t10 - 2 .* t11 + t5) - 2 .* t29 .* t28 .* x_m(vecBI)) + t9 .* (t16 - t18 + t10 + t11) + 2 .* t7 .* t36 + (t6 + t11) .* t58) .* (t23 .* ((xi_3(vecBI) - x_m(vecBI)) .* xi_2(vecBI) - t9 + t7 - t29 .* t28) + (-t21 + (xi_3(vecBI) + x_m(vecBI)) .* xi_2(vecBI) - t7 + t36) .* t13) .* z_m(vecBI) ./ t23 ./ t13 .* t2;
    t1 = (eta_2(vecBI) .^ 2);
    t2 = (y_m(vecBI) .* eta_2(vecBI));
    t4 = (x_m(vecBI) .^ 2);
    t7 = (xi_2(vecBI) .^ 2);
    t8 = (y_m(vecBI) .^ 2);
    t9 = (z_m(vecBI) .^ 2);
    t11 = sqrt((-2 .* x_m(vecBI) .* xi_2(vecBI) + t1 - 2 .* t2 + t4 + t7 + t8 + t9));
    t13 = (eta_3(vecBI) .^ 2);
    t14 = (y_m(vecBI) .* eta_3(vecBI));
    t16 = x_m(vecBI) .* xi_3(vecBI);
    t18 = xi_3(vecBI) .^ 2;
    t20 = sqrt((t13 - 2 .* t14 + t4 - 2 .* t16 + t18 + t8 + t9));
    t23 = (C(vecBI) .^ 2);
    t34 = (D_LE(vecBI) .^ 2);
    t38 = (-y_m(vecBI) + eta_3(vecBI));
    t40 = (y_m(vecBI) - eta_2(vecBI));
    t42 = (eta_2(vecBI) - eta_3(vecBI));
    t43 = t42 .* x_m(vecBI);
    t47 = -t38;
    t53 = t40 .* t42;
    t60 = 12 .* t14;
    t61 = 6 .* t8;
    t62 = 6 .* t9;
    t81 = t42 .^ 2;
    tmp31 = -0.1e1 ./ (t7 .* (6 .* t13 - t60 + t61 + t62) + xi_2(vecBI) .* (xi_3(vecBI) .* (12 .* eta_2(vecBI) .* t47 + t60 - 12 .* t8 - 12 .* t9) - 12 .* t47 .* t43) + t18 .* (6 .* t1 - 12 .* t2 + t61 + t62) + 12 .* t16 .* t53 + 6 .* (t4 + t9) .* t81) .* (-xi_3(vecBI) + xi_2(vecBI)) .* (t11 .* ((x_m(vecBI) - xi_3(vecBI)) .* xi_2(vecBI) + t18 - t16 + t47 .* t42) - (-t7 + (xi_3(vecBI) + x_m(vecBI)) .* xi_2(vecBI) - t16 + t53) .* t20) .* (xi_2(vecBI) .* t38 + xi_3(vecBI) .* t40 + t43) .* (t7 .* t23 + C(vecBI) .* (C(vecBI) .* xi_3(vecBI) + 3 .* D_LE(vecBI)) .* xi_2(vecBI) + t18 .* t23 + 3 .* C(vecBI) .* D_LE(vecBI) .* xi_3(vecBI) + 3 .* t34) ./ t20 ./ t11;
    t1 = (eta_2(vecBI) .^ 2);
    t2 = (y_m(vecBI) .* eta_2(vecBI));
    t4 = (x_m(vecBI) .^ 2);
    t7 = (xi_2(vecBI) .^ 2);
    t8 = (y_m(vecBI) .^ 2);
    t9 = (z_m(vecBI) .^ 2);
    t11 = sqrt((-2 .* x_m(vecBI) .* xi_2(vecBI) + t1 - 2 .* t2 + t4 + t7 + t8 + t9));
    t13 = (eta_3(vecBI) .^ 2);
    t14 = (y_m(vecBI) .* eta_3(vecBI));
    t16 = x_m(vecBI) .* xi_3(vecBI);
    t18 = xi_3(vecBI) .^ 2;
    t20 = sqrt((t13 - 2 .* t14 + t4 - 2 .* t16 + t18 + t8 + t9));
    t23 = (-y_m(vecBI) + eta_3(vecBI));
    t25 = (y_m(vecBI) - eta_2(vecBI));
    t27 = (eta_2(vecBI) - eta_3(vecBI));
    t28 = t27 .* x_m(vecBI);
    t33 = -t23;
    t39 = t25 .* t27;
    t50 = 4 .* t14;
    t51 = 2 .* t8;
    t52 = 2 .* t9;
    t71 = t27 .^ 2;
    tmp32 = -0.1e1 ./ (t7 .* (2 .* t13 - t50 + t51 + t52) + xi_2(vecBI) .* (xi_3(vecBI) .* (4 .* eta_2(vecBI) .* t33 + t50 - 4 .* t8 - 4 .* t9) - 4 .* t33 .* t28) + t18 .* (2 .* t1 - 4 .* t2 + t51 + t52) + 4 .* t16 .* t39 + 2 .* (t4 + t9) .* t71) .* (-xi_3(vecBI) + xi_2(vecBI)) .* (C(vecBI) .* xi_2(vecBI) + C(vecBI) .* xi_3(vecBI) + 2 .* D_LE(vecBI)) .* (t11 .* ((x_m(vecBI) - xi_3(vecBI)) .* xi_2(vecBI) + t18 - t16 + t33 .* t27) - (-t7 + (xi_3(vecBI) + x_m(vecBI)) .* xi_2(vecBI) - t16 + t39) .* t20) .* (xi_2(vecBI) .* t23 + xi_3(vecBI) .* t25 + t28) ./ t20 ./ t11;
    t2 = (eta_3(vecBI) .^ 2);
    t3 = (y_m(vecBI) .* eta_3(vecBI));
    t5 = (x_m(vecBI) .^ 2);
    t6 = (x_m(vecBI) .* xi_3(vecBI));
    t8 = (xi_3(vecBI) .^ 2);
    t9 = (y_m(vecBI) .^ 2);
    t10 = (z_m(vecBI) .^ 2);
    t12 = sqrt((t2 - 2 .* t3 + t5 - 2 .* t6 + t8 + t9 + t10));
    t15 = (-y_m(vecBI) + eta_3(vecBI));
    t17 = (y_m(vecBI) - eta_2(vecBI));
    t19 = (eta_2(vecBI) - eta_3(vecBI));
    t20 = (t19 .* x_m(vecBI));
    t23 = (eta_2(vecBI) .^ 2);
    t24 = (y_m(vecBI) .* eta_2(vecBI));
    t28 = (xi_2(vecBI) .^ 2);
    t30 = sqrt((-2 .* x_m(vecBI) .* xi_2(vecBI) + t10 + t23 - 2 .* t24 + t28 + t5 + t9));
    t34 = -t15;
    t40 = (t17 .* t19);
    t48 = 12 .* t3;
    t49 = 6 .* t9;
    t50 = 6 .* t10;
    t69 = (t19 .^ 2);
    tmp33 = 0.1e1 ./ (t28 .* (6 .* t2 - t48 + t49 + t50) + xi_2(vecBI) .* (xi_3(vecBI) .* (12 .* eta_2(vecBI) .* t34 - 12 .* t10 + t48 - 12 .* t9) - 12 .* t34 .* t20) + t8 .* (6 .* t23 - 12 .* t24 + t49 + t50) + 12 .* t6 .* t40 + 6 .* (t5 + t10) .* t69) .* (xi_3(vecBI) .* xi_2(vecBI) + t28 + t8) .* (t30 .* ((xi_3(vecBI) - x_m(vecBI)) .* xi_2(vecBI) - t8 + t6 - t34 .* t19) + (-t28 + (xi_3(vecBI) + x_m(vecBI)) .* xi_2(vecBI) - t6 + t40) .* t12) ./ t30 .* (xi_2(vecBI) .* t15 + xi_3(vecBI) .* t17 + t20) ./ t12 .* (-xi_3(vecBI) + xi_2(vecBI));
    t2 = (eta_3(vecBI) .^ 2);
    t3 = (y_m(vecBI) .* eta_3(vecBI));
    t5 = (x_m(vecBI) .^ 2);
    t6 = (x_m(vecBI) .* xi_3(vecBI));
    t8 = (xi_3(vecBI) .^ 2);
    t9 = (y_m(vecBI) .^ 2);
    t10 = (z_m(vecBI) .^ 2);
    t12 = sqrt((t2 - 2 .* t3 + t5 - 2 .* t6 + t8 + t9 + t10));
    t15 = (-y_m(vecBI) + eta_3(vecBI));
    t17 = (y_m(vecBI) - eta_2(vecBI));
    t19 = (eta_2(vecBI) - eta_3(vecBI));
    t20 = (t19 .* x_m(vecBI));
    t24 = (eta_2(vecBI) .^ 2);
    t25 = (y_m(vecBI) .* eta_2(vecBI));
    t29 = (xi_2(vecBI) .^ 2);
    t31 = sqrt((-2 .* x_m(vecBI) .* xi_2(vecBI) + t10 + t24 - 2 .* t25 + t29 + t5 + t9));
    t36 = -t15;
    t42 = (t17 .* t19);
    t47 = 4 .* t3;
    t48 = 2 .* t9;
    t49 = 2 .* t10;
    t68 = (t19 .^ 2);
    tmp34 = 0.1e1 ./ (t29 .* (2 .* t2 - t47 + t48 + t49) + xi_2(vecBI) .* (xi_3(vecBI) .* (4 .* eta_2(vecBI) .* t36 - 4 .* t10 + t47 - 4 .* t9) - 4 .* t36 .* t20) + t8 .* (2 .* t24 - 4 .* t25 + t48 + t49) + 4 .* t6 .* t42 + 2 .* (t5 + t10) .* t68) .* (t31 .* ((xi_3(vecBI) - x_m(vecBI)) .* xi_2(vecBI) - t8 + t6 - t36 .* t19) + (-t29 + (xi_3(vecBI) + x_m(vecBI)) .* xi_2(vecBI) - t6 + t42) .* t12) ./ t31 .* (xi_2(vecBI) + xi_3(vecBI)) .* (xi_2(vecBI) .* t15 + xi_3(vecBI) .* t17 + t20) ./ t12 .* (-xi_3(vecBI) + xi_2(vecBI));
    t1 = (eta_2(vecBI) .^ 2);
    t2 = (y_m(vecBI) .* eta_2(vecBI));
    t4 = (x_m(vecBI) .^ 2);
    t7 = (xi_2(vecBI) .^ 2);
    t8 = (y_m(vecBI) .^ 2);
    t9 = (z_m(vecBI) .^ 2);
    t11 = sqrt((-2 .* x_m(vecBI) .* xi_2(vecBI) + t1 - 2 .* t2 + t4 + t7 + t8 + t9));
    t13 = (eta_3(vecBI) .^ 2);
    t14 = (y_m(vecBI) .* eta_3(vecBI));
    t16 = x_m(vecBI) .* xi_3(vecBI);
    t18 = xi_3(vecBI) .^ 2;
    t20 = sqrt((t13 - 2 .* t14 + t4 - 2 .* t16 + t18 + t8 + t9));
    t23 = (-y_m(vecBI) + eta_3(vecBI));
    t25 = (y_m(vecBI) - eta_2(vecBI));
    t27 = (eta_2(vecBI) - eta_3(vecBI));
    t28 = t27 .* x_m(vecBI);
    t33 = -t23;
    t39 = t25 .* t27;
    t48 = (C(vecBI) .* xi_3(vecBI)) + 0.3e1 ./ 0.2e1 .* D_LE(vecBI);
    t53 = 12 .* t14;
    t54 = 6 .* t8;
    t55 = 6 .* t9;
    t74 = t27 .^ 2;
    tmp35 = -0.2e1 ./ (t7 .* (6 .* t13 - t53 + t54 + t55) + xi_2(vecBI) .* (xi_3(vecBI) .* (12 .* eta_2(vecBI) .* t33 + t53 - 12 .* t8 - 12 .* t9) - 12 .* t33 .* t28) + t18 .* (6 .* t1 - 12 .* t2 + t54 + t55) + 12 .* t16 .* t39 + 6 .* (t4 + t9) .* t74) .* ((t7 .* C(vecBI)) + xi_2(vecBI) .* t48 + t48 .* xi_3(vecBI)) .* (-xi_3(vecBI) + xi_2(vecBI)) .* (t11 .* ((x_m(vecBI) - xi_3(vecBI)) .* xi_2(vecBI) + t18 - t16 + t33 .* t27) - (-t7 + (xi_3(vecBI) + x_m(vecBI)) .* xi_2(vecBI) - t16 + t39) .* t20) .* (xi_2(vecBI) .* t23 + xi_3(vecBI) .* t25 + t28) ./ t20 ./ t11;
    t2 = (eta_3(vecBI) .^ 2);
    t4 = 2 .* y_m(vecBI) .* eta_3(vecBI);
    t5 = (x_m(vecBI) .^ 2);
    t6 = (x_m(vecBI) .* xi_3(vecBI));
    t8 = (xi_3(vecBI) .^ 2);
    t9 = y_m(vecBI) .^ 2;
    t10 = (z_m(vecBI) .^ 2);
    t12 = sqrt((t2 - t4 + t5 - 2 .* t6 + t8 + t9 + t10));
    t15 = -y_m(vecBI) + eta_3(vecBI);
    t17 = y_m(vecBI) - eta_2(vecBI);
    t19 = eta_2(vecBI) - eta_3(vecBI);
    t20 = t19 .* x_m(vecBI);
    t23 = eta_2(vecBI) .^ 2;
    t25 = 2 .* y_m(vecBI) .* eta_2(vecBI);
    t28 = (xi_2(vecBI) .^ 2);
    t30 = sqrt((-2 .* x_m(vecBI) .* xi_2(vecBI) + t10 + t23 - t25 + t28 + t5 + t9));
    t34 = -t15;
    t40 = t17 .* t19;
    t61 = t19 .^ 2;
    tmp36 = 0.1e1 ./ (t28 .* (t2 - t4 + t9 + t10) + xi_2(vecBI) .* (xi_3(vecBI) .* (2 .* eta_2(vecBI) .* t34 - 2 .* t10 + t4 - 2 .* t9) - 2 .* t34 .* t20) + t8 .* (t23 - t25 + t9 + t10) + 2 .* t6 .* t40 + (t5 + t10) .* t61) .* (t30 .* ((xi_3(vecBI) - x_m(vecBI)) .* xi_2(vecBI) - t8 + t6 - t34 .* t19) + (-t28 + (xi_3(vecBI) + x_m(vecBI)) .* xi_2(vecBI) - t6 + t40) .* t12) ./ t30 .* (xi_2(vecBI) .* t15 + xi_3(vecBI) .* t17 + t20) ./ t12 .* (-xi_3(vecBI) + xi_2(vecBI));
    
    eps = 1e-5;
    r1 = [xi_2(vecBI) eta_2(vecBI) xi_2(vecBI).*0] - [x_m(vecBI) y_m(vecBI) z_m(vecBI)];
    r2 = [xi_3(vecBI) eta_3(vecBI) xi_3(vecBI).*0] - [x_m(vecBI) y_m(vecBI) z_m(vecBI)];
    idx_corr = sqrt(sum(r1.^2,2)) < eps | sqrt(sum(r2.^2,2)) < eps | sqrt(sum(cross(r1,r2,2).^2,2)) < eps;
    
    tmp11(idx_corr) = tmp11(idx_corr).*0;
    tmp12(idx_corr) = tmp11(idx_corr).*0;
    tmp13(idx_corr) = tmp11(idx_corr).*0;
    tmp14(idx_corr) = tmp11(idx_corr).*0;
    tmp15(idx_corr) = tmp11(idx_corr).*0;
    tmp16(idx_corr) = tmp11(idx_corr).*0;
    tmp21(idx_corr) = tmp11(idx_corr).*0;
    tmp22(idx_corr) = tmp11(idx_corr).*0;
    tmp23(idx_corr) = tmp11(idx_corr).*0;
    tmp24(idx_corr) = tmp11(idx_corr).*0;
    tmp25(idx_corr) = tmp11(idx_corr).*0;
    tmp26(idx_corr) = tmp11(idx_corr).*0;
    tmp31(idx_corr) = tmp11(idx_corr).*0;
    tmp32(idx_corr) = tmp11(idx_corr).*0;
    tmp33(idx_corr) = tmp11(idx_corr).*0;
    tmp34(idx_corr) = tmp11(idx_corr).*0;
    tmp35(idx_corr) = tmp11(idx_corr).*0;
    tmp36(idx_corr) = tmp11(idx_corr).*0;
    
end

%%
infl_new(1,1,:) = reshape(tmp11,1,1,[]);
infl_new(1,2,:) = reshape(tmp12,1,1,[]);
infl_new(1,3,:) = reshape(tmp13,1,1,[]);
infl_new(1,4,:) = reshape(tmp14,1,1,[]);
infl_new(1,5,:) = reshape(tmp15,1,1,[]);
infl_new(1,6,:) = reshape(tmp16,1,1,[]);

infl_new(2,1,:) = reshape(tmp21,1,1,[]);
infl_new(2,2,:) = reshape(tmp22,1,1,[]);
infl_new(2,3,:) = reshape(tmp23,1,1,[]);
infl_new(2,4,:) = reshape(tmp24,1,1,[]);
infl_new(2,5,:) = reshape(tmp25,1,1,[]);
infl_new(2,6,:) = reshape(tmp26,1,1,[]);

infl_new(3,1,:) = reshape(tmp31,1,1,[]);
infl_new(3,2,:) = reshape(tmp32,1,1,[]);
infl_new(3,3,:) = reshape(tmp33,1,1,[]);
infl_new(3,4,:) = reshape(tmp34,1,1,[]);
infl_new(3,5,:) = reshape(tmp35,1,1,[]);
infl_new(3,6,:) = reshape(tmp36,1,1,[]);

%%
infl_loc = real(infl_new);
infl_loc(:,:,idx_flp) = -infl_loc(:,:,idx_flp);

idx_nan = find(reshape(sum(any(isnan(infl_new) | isinf(infl_new))),[],1) > 0);
if any(idx_nan)
    %     disp('Nan induction in fcnHDVEIND_DB');
end
infl_loc(isnan(infl_loc) | isinf(infl_loc)) = 0;

% infl_loc(1:2,:,idx_afx) = infl_loc(1:2,:,idx_afx).*reshape((z_m_orig./ztol),1,1,[]);

end