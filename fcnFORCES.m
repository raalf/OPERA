function [CL, CDi, CY, e, vecDVELIFT, vecDVEDRAG] = fcnFORCES(valTIMESTEP, matVLST, matCENTER, matELST, matROTANG, matUINF, matCOEFF, vecTEDVE, valDENSITY, valNELE, matSPANDIR, vecTE, vecDVEAREA, matPLEX)

vecDVELIFT = nan(valNELE,1);
vecDVEDRAG = nan(valNELE,1);
vecDVESIDE = nan(valNELE,1);
matDVEDRAG_DIR = matUINF./sqrt(sum(matUINF.^2,2));
matDVELIFT_DIR = cross(matDVEDRAG_DIR, matSPANDIR, 2);
matDVESIDE_DIR = cross(matDVELIFT_DIR, matDVEDRAG_DIR, 2);

liftfree = nan(valNELE,1);
sidefree = nan(valNELE,1);
liftind = nan(valNELE,1);
sideind = nan(valNELE,1);
dragind = nan(valNELE,1);

%% To element local
xi_1 = permute(matPLEX(1,1,vecTEDVE),[3 2 1]);
xi_2 = permute(matPLEX(2,1,vecTEDVE),[3 2 1]);
xi_3 = permute(matPLEX(3,1,vecTEDVE),[3 2 1]);

eta_1 = permute(matPLEX(1,2,vecTEDVE),[3 2 1]);
eta_2 = permute(matPLEX(2,2,vecTEDVE),[3 2 1]);
eta_3 = permute(matPLEX(3,2,vecTEDVE),[3 2 1]);

idx_flp = abs(xi_2 - xi_3) < 1e-3 & abs(xi_1 - xi_2) > 1e-3;
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

%% Freestream Forces (via Kutty J)
% Integrating the cross product of the circulation along the TE and the
% density times U_INF. Circulation equation is put solely in terms of xi
% by defining eta in terms of xi. It is integrated from xi_1 to xi_3.
% T.D.K 2019-01-23 230 KING ST. EAST, TORONTO, ON, CANADA
uvw = fcnGLOBSTAR(matUINF(vecTEDVE,:), matROTANG(vecTEDVE,:));

% For readability:
A_1 = matCOEFF(vecTEDVE,1); A_2 = matCOEFF(vecTEDVE,2); B_1 = matCOEFF(vecTEDVE,3);
B_2 = matCOEFF(vecTEDVE,4); C_3 = matCOEFF(vecTEDVE,5); 
u = uvw(:,1); v = uvw(:,2); w = uvw(:,3);

% Analytically integrating circulation along the TE vector
% Projected circulation onto trailing edge
% term = ((A_1.*E.^2 + B_1).*xi_1.^2 + ((A_1.*E.^2 + B_1).*xi_3 + (3.*A_1.*D_TE + 3.*A_2).*E + 3.*B_2).*xi_1 + (A_1.*E.^2 + B_1).*xi_3.^2 + ((3.*A_1.*D_TE + 3.*A_2).*E + 3.*B_2).*xi_3 + 3.*A_1.*D_TE.^2 + 6.*A_2.*D_TE + 6.*C_3);
% Spanwise circulation
term = (B_1.*xi_1.^2 + (B_1.*xi_3 + 3.*B_2).*xi_1 + B_1.*xi_3.^2 + 3.*B_2.*xi_3 + 6.*C_3);
tmp21 = -(1/6).*((xi_1 - xi_3).*valDENSITY.*term.*w.*(eta_1 - eta_3));
tmp22 =  (1/6).*((xi_1 - xi_3).^2.*valDENSITY.*term.*w.*(eta_1 - eta_3));
tmp33 = -(1/6).*(xi_1.*v - xi_3.*v - u.*(eta_3 - eta_1)).*term.*(xi_1 - xi_3).*valDENSITY;

% This is the special force
F = fcnSTARGLOB(-[tmp21 tmp22 tmp33]./(xi_3 - xi_1), matROTANG(vecTEDVE,:));

% Splitting special force into lift and side forces
liftfree(vecTEDVE,1) = dot(F, matDVELIFT_DIR(vecTEDVE,:), 2);
sidefree(vecTEDVE,1) = dot(F, matDVESIDE_DIR(vecTEDVE,:), 2);

vecDVELIFT(vecTEDVE,1) = liftfree(vecTEDVE,1);
vecDVESIDE(vecTEDVE,1) = sidefree(vecTEDVE,1);

%% Induced Forces 


%% Output
CL = nansum(vecDVELIFT)./(0.5.*valDENSITY.*sum(vecDVEAREA));
CDi = nansum(vecDVEDRAG)./(0.5.*valDENSITY.*sum(vecDVEAREA));
CY = nansum(vecDVESIDE)./(0.5.*valDENSITY.*sum(vecDVEAREA));
e = (CL.^2)./(pi.*((max(matVLST(:,2))-min(matVLST(:,2))).^2)./sum(vecDVEAREA).*CDi);
fprintf('Timestep: %d\t\tCL = %0.5f\t\tCDi = %0.5f\t\te = %0.5f\n', valTIMESTEP, CL, CDi, e);

end