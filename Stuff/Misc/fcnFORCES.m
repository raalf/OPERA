function [CL, CDi, CY, e, vecDVELIFT, vecDVEDRAG, matDVEDRAG_DIR, matDVELIFT_DIR, matDVESIDE_DIR] = fcnWFORCES(valTIMESTEP, matVLST, matCENTER, matELST, matROTANG, matUINF, matCOEFF, vecTEDVE, valDENSITY, valNELE, matSPANDIR, vecTE, vecDVEAREA, matPLEX, matWCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matVUINF, matWVLST, vecWLE, vecWLEDVE, matWELST, valAREA, valSPAN, matWDVECT, matDVECT, vecDVESYM, vecWDVESYM)
lim = 1e10;

%% Initializing
vecDVELIFT = nan(valNELE,1);
vecDVEDRAG = nan(valNELE,1);
vecDVESIDE = nan(valNELE,1);
% matDVEDRAG_DIR = matUINF./sqrt(sum(matUINF.^2,2));
% matDVELIFT_DIR = cross(matDVEDRAG_DIR, matSPANDIR, 2);
% matDVESIDE_DIR = cross(matDVELIFT_DIR, matDVEDRAG_DIR, 2);


matDVEDRAG_DIR = repmat(matUINF(1,:)./sqrt(sum(matUINF(1,:).^2,2)), valNELE, 1);
matDVELIFT_DIR = repmat(cross(matDVEDRAG_DIR(1,:), matSPANDIR(1,:), 2), valNELE, 1);
matDVESIDE_DIR = repmat(cross(matDVELIFT_DIR(1,:), matDVEDRAG_DIR(1,:), 2), valNELE,1);

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

C = (eta_3 - eta_2)./(xi_3 - xi_2);
D_LE = eta_2 - ((xi_2.*(eta_3 - eta_2))./(xi_3 - xi_2));
E = (eta_3 - eta_1)./(xi_3 - xi_1);
D_TE = eta_1 - ((xi_1.*(eta_3 - eta_1))./(xi_3 - xi_1));

%% Freestream Forces (via Kutty J)
% Integrating the cross product of the circulation along the TE and U_INF.
% T.D.K 2019-01-23 230 KING ST. EAST, TORONTO, ON, CANADA
% tescale = sqrt(sum([xi_3 - xi_1 eta_3 - eta_1].^2,2))./(xi_3 - xi_1);

% Velocity along the TE
locations = [0.2:0.3:0.8]';
uvw_te = locations.*reshape(matVUINF(matELST(vecTE,1),:)', 1, 3, []) + (1 - locations).*reshape(matVUINF(matELST(vecTE,2),:)', 1, 3, []);
fpl_te = locations.*reshape([xi_1 eta_1 xi_1.*0]', 1, 3, []) + (1 - locations).*reshape([xi_3 eta_3 xi_3.*0]', 1, 3, []);

tau = fpl_te(:,1,:);
for i = 1:size(vecTEDVE,1)
    uvw = fcnGLOBSTAR(uvw_te(:,:,i), repmat(matROTANG(vecTEDVE(i),:), length(locations), 1));
    u(i,:) = polyfit(tau(:,1,i), uvw(:,1), 2);
    v(i,:) = polyfit(tau(:,1,i), uvw(:,2), 2);
    w(i,:) = polyfit(tau(:,1,i), uvw(:,3), 2);    
end

A_1 = matCOEFF(vecTEDVE,1); A_2 = matCOEFF(vecTEDVE,2); B_1 = matCOEFF(vecTEDVE,3);
B_2 = matCOEFF(vecTEDVE,4); C_2 = matCOEFF(vecTEDVE,5); C_3 = matCOEFF(vecTEDVE,6);

u_out = fcnKJU(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY, lim);
v_out = fcnKJV(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY, lim);
w_out = fcnKJW(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY, lim);
F = [u_out v_out w_out];

F = fcnSTARGLOB(-F, matROTANG(vecTEDVE,:));

% Splitting special force into lift and side forces
liftfree(vecTEDVE,1) = dot(F, matDVELIFT_DIR(vecTEDVE,:), 2);
sidefree(vecTEDVE,1) = dot(F, matDVESIDE_DIR(vecTEDVE,:), 2);

%% Induced Forces
% Induced velocities at wake leading edge DVEs (wind is 3 x 3 x num_dve) 
fpg_og = reshape(locations,1,1,[]).*matWVLST(matWELST(vecWLE,1),:) + (1 - reshape(locations,1,1,[])).*matWVLST(matWELST(vecWLE,2),:);
fpg_og = reshape(permute(fpg_og, [2 1 3]), size(fpg_og, 2), [])';
vecBOUNDIND = false(valWNELE,1);
% vecBOUNDIND(vecWLEDVE) = true;
wind = fcnSDVEVEL(fpg_og, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, vecWDVESYM, vecBOUNDIND, 0);
% hold on
% quiver3(fpg_og(:,1), fpg_og(:,2), fpg_og(:,3), wind(:,1), wind(:,2), wind(:,3), 'b');
% hold off
wind = permute(reshape(wind', 3, [], 3), [3 1 2]);

% Velocity along the TE
u = []; v = []; w = [];
for i = 1:size(vecTEDVE,1)
    uvw = fcnGLOBSTAR(wind(:,:,i), repmat(matROTANG(vecTEDVE(i),:), length(locations), 1));
    u(i,:) = polyfit(tau(:,1,i), uvw(:,1), 2);
    v(i,:) = polyfit(tau(:,1,i), uvw(:,2), 2);
    w(i,:) = polyfit(tau(:,1,i), uvw(:,3), 2);   
end

u_out = fcnKJU(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY, lim);
v_out = fcnKJV(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY, lim);
w_out = fcnKJW(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY, lim);
F = [u_out v_out w_out];

F = fcnSTARGLOB(-F, matROTANG(vecTEDVE,:));

% Splitting induced force into lift and side forces
liftind(vecTEDVE,1) = dot(F, matDVELIFT_DIR(vecTEDVE,:), 2);
sideind(vecTEDVE,1) = dot(F, matDVESIDE_DIR(vecTEDVE,:), 2);
dragind(vecTEDVE,1) = dot(F, matDVEDRAG_DIR(vecTEDVE,:), 2);

% vecDVELIFT = liftfree + liftind;
vecDVELIFT = liftfree;
vecDVESIDE = sidefree + sideind;
vecDVEDRAG = dragind;

% Symmetry
vecDVELIFT(vecDVESYM) = vecDVELIFT(vecDVESYM)*2;
vecDVEDRAG(vecDVESYM) = vecDVEDRAG(vecDVESYM)*2;

%% Output
CL = nansum(vecDVELIFT)./(0.5.*valDENSITY.*valAREA);
CDi = nansum(vecDVEDRAG)./(0.5.*valDENSITY.*valAREA);
CY = nansum(vecDVESIDE)./(0.5.*valDENSITY.*valAREA);
% e = (CL.^2)./(pi.*((valSPAN.^2)./valAREA).*CDi);
q_inf = 0.5.*valDENSITY;
e = ((nansum(vecDVELIFT)./q_inf).^2)./(pi.*(valSPAN.^2).*(nansum(vecDVEDRAG)./q_inf));

fprintf('Timestep: %d\t\tCL = %0.5f\t\tCDi = %0.5f\t\te = %0.5f\n', valTIMESTEP, CL, CDi, e);

end