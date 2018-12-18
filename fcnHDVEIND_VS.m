function [infl_loc] = fcnHDVEIND_VS(dvenum, dvetype, fpg, matPLEX, matROTANG, matCONTROL)
warning('on')
tol = 1e-10;
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

idx_flp = abs(xi_2 - xi_3) < tol;
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

margin_edge = 1e-5;
margin_on_element = 1e-5;
xi_left = min([xi_1, xi_3],[],2);
xi_right = max([xi_1, xi_3],[],2);
idx_on_element = y_m >= te_eta - margin_edge & y_m <= le_eta + margin_edge & x_m >= xi_left - margin_edge & x_m <= xi_right + margin_edge & abs(z_m) <= margin_on_element;
idx_on_edge =   (abs(y_m - te_eta) < margin_edge & (xi_left - margin_edge <= x_m & x_m <= xi_right + margin_edge) & abs(z_m) <= margin_on_element) | ...
                (abs(y_m - le_eta) < margin_edge & (xi_left - margin_edge <= x_m & x_m <= xi_right + margin_edge) & abs(z_m) <= margin_on_element) | ...
                (abs(x_m - xi_left) < margin_edge & (te_eta - margin_edge <= y_m & y_m <= le_eta + margin_edge) & abs(z_m) <= margin_on_element) | ...
                (abs(x_m - xi_right) < margin_edge & (te_eta - margin_edge <= y_m & y_m <= le_eta + margin_edge) & abs(z_m) <= margin_on_element);
            
% disp(['Edge calls: ', num2str(sum(idx_on_edge))]);
         
%%
infl_new = zeros(3,5,len);
k = repmat(0.1, len, 1);
hchord = ones(len, 1);

%% % 2-3 minus 1-3
% 2-3
hspan = abs(xi_3 - xi_2)./2;
phi = atan(-C);
xsiA = fpl - [(xi_2 + xi_3)./2 (eta_2 + eta_3)./2 zeros(len,1)];
xsiA = fcnGLOBSTAR(xsiA, repmat([0 0 pi/2], len, 1));
[~, btmp, ctmp] = fcnVSIND(hspan, hchord, phi, xsiA, k, 0);
bloc23 = fcnSTARGLOB(btmp, repmat([0 0 pi/2],len,1));
cloc23 = fcnSTARGLOB(ctmp, repmat([0 0 pi/2],len,1));

% 1-3
hspan = abs(xi_3 - xi_1)./2;
phi = atan(-E);
xsiA = fpl - [(xi_1 + xi_3)./2 (eta_1 + eta_3)./2 zeros(len,1)];
xsiA = fcnGLOBSTAR(xsiA, repmat([0 0 pi/2], len, 1));
[~, btmp, ctmp] = fcnVSIND(hspan, hchord, phi, xsiA, k, 0);
bloc13 = fcnSTARGLOB(btmp, repmat([0 0 pi/2],len,1));
cloc13 = fcnSTARGLOB(ctmp, repmat([0 0 pi/2],len,1));
infl_new(:,3,:) = reshape([cloc23-cloc13]',3,size(len,1),[]);
infl_new(:,4,:) = reshape([bloc23-bloc13]',3,size(len,1),[]);

%% 
% 1-2
hspan = abs(eta_2 - eta_1)./2;
phi = zeros(len,1);
xsiA = fpl - [(xi_1 + xi_2)./2 (eta_1 + eta_2)./2 zeros(len,1)];
[~, bloc12, cloc12] = fcnVSIND(hspan, hchord, phi, xsiA, k, 0);

% 2-3
hspan = abs(eta_3 - eta_2)./2;
phi = atan(1./C);
xsiA = fpl - [(xi_2 + xi_3)./2 (eta_2 + eta_3)./2 zeros(len,1)];
[~, bloc23, cloc23] = fcnVSIND(hspan, hchord, phi, xsiA, k, 0);

% 1-3
hspan = abs(eta_3 - eta_1)./2;
phi = atan(1./E);
xsiA = fpl - [(xi_1 + xi_3)./2 (eta_1 + eta_3)./2 zeros(len,1)];
[~, bloc13, cloc13] = fcnVSIND(hspan, hchord, phi, xsiA, k, 0);

% Compiling
idx = eta_3 > eta_2 + tol;
infl_new(:,1,idx) = reshape([cloc12(idx,:) + cloc23(idx,:) - cloc13(idx,:)]',3,size(len,1),[]);
infl_new(:,2,idx) = reshape([bloc12(idx,:) + bloc23(idx,:) - bloc13(idx,:)]',3,size(len,1),[]);

idx = abs(eta_3 - eta_2) <= tol;
infl_new(:,1,idx) = reshape([cloc12(idx,:) - cloc13(idx,:)]',3,size(len,1),[]);
infl_new(:,2,idx) = reshape([bloc12(idx,:) - bloc13(idx,:)]',3,size(len,1),[]);

idx = (eta_3 < eta_2 - tol) & (eta_3 > eta_1 + tol);
infl_new(:,1,idx) = reshape([cloc12(idx,:) - cloc23(idx,:) - cloc13(idx,:)]',3,size(len,1),[]);
infl_new(:,2,idx) = reshape([bloc12(idx,:) - bloc23(idx,:) - bloc13(idx,:)]',3,size(len,1),[]);

idx = abs(eta_3 - eta_1) <= tol;
infl_new(:,1,idx) = reshape([cloc12(idx,:) - cloc23(idx,:)]',3,size(len,1),[]);
infl_new(:,2,idx) = reshape([bloc12(idx,:) - bloc23(idx,:)]',3,size(len,1),[]);

idx = eta_3 < eta_1 - tol;
infl_new(:,1,idx) = reshape([cloc12(idx,:) - cloc23(idx,:) + cloc13(idx,:)]',3,size(len,1),[]);
infl_new(:,2,idx) = reshape([bloc12(idx,:) - bloc23(idx,:) + bloc13(idx,:)]',3,size(len,1),[]);

%%
infl_loc = real(infl_new);
infl_loc(:,:,idx_flp) = -infl_loc(:,:,idx_flp);

idx_nan = find(reshape(sum(any(isnan(infl_new) | isinf(infl_new))),[],1) > 0);
% disp(['Inf or NaN induction: ', num2str(length(idx_nan))]);

% infl_loc(isnan(infl_loc) | isinf(infl_loc)) = 0;
end