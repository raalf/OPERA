function [vecDVELIFT, vecDVEDRAG, vecDVESIDE, matDVEDRAG_DIR, matDVELIFT_DIR, matDVESIDE_DIR, vecDGAMMA_DT, vecDGAMMA_DETA] = fcnDVEFORCES(valTIMESTEP, strWAKE_TYPE, matVLST, matELST, matROTANG, ...
    matUINF, matCOEFF, vecTEDVE, valDENSITY, valNELE, matSPANDIR, vecTE, matPLEX, matWCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG, ...
    matVUINF, matWVLST, vecWLE, vecWLEDVE, matWELST, valAREA, valSPAN, vecDVESYM, vecWDVESYM, vecDGAMMA_DT, vecDGAMMA_DETA, valWSIZE, matWDVE, vecWDVEFLIP, matWEIDX, vecWEMU, vecWVMU)
lim = 1e10;

%% Initializing
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

C = (eta_3 - eta_2)./(xi_3 - xi_2);
D_LE = eta_2 - ((xi_2.*(eta_3 - eta_2))./(xi_3 - xi_2));
E = (eta_3 - eta_1)./(xi_3 - xi_1);
D_TE = eta_1 - ((xi_1.*(eta_3 - eta_1))./(xi_3 - xi_1));

%% Freestream Forces (via Kutty J)
% Integrating the cross product of the circulation along the TE and U_INF.
% T.D.K 2019-01-23 230 KING ST. EAST, TORONTO, ON, CANADA

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

% Freestream lift and side force (U x spanwise_direction)
tmp = fcnGLOBSTAR(matSPANDIR(vecTEDVE,:), matROTANG(vecTEDVE,:));
u_out = fcnKJU(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY, tmp(:,1), tmp(:,2), tmp(:,3), lim);
v_out = fcnKJV(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY, tmp(:,1), tmp(:,2), tmp(:,3), lim);
w_out = fcnKJW(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY, tmp(:,1), tmp(:,2), tmp(:,3), lim);
F = fcnSTARGLOB([u_out v_out w_out], matROTANG(vecTEDVE,:));

% Splitting special force into lift and side forces
liftfree(vecTEDVE,1) = dot(F, matDVELIFT_DIR(vecTEDVE,:), 2);
sidefree(vecTEDVE,1) = dot(F, matDVESIDE_DIR(vecTEDVE,:), 2);

%% Induced Forces
% Induced velocities at wake leading edge DVEs (wind is 3 x 3 x num_dve)

locations = linspace(0.2, 0.8, 3)';

wind = nan(size(locations,1), 3, size(vecWLEDVE,1));
fpg_og = nan(size(locations,1), 3, size(vecWLEDVE,1));

boundind = false(valWNELE,1);
boundind(vecWLEDVE) = true;

% Need to adjust wake LE elements to straight, one by one
for i = 1:valWSIZE
    tmpWVLST = matWVLST;
    
    tmp_verts(1) = matWDVE(vecWLEDVE(i),2);
    tmp_verts(2) = matWDVE(vecWLEDVE(i),3);
    
    tmp_mp = mean(matWVLST(tmp_verts,:),1);
    tmp_spandir = matSPANDIR(vecTEDVE(i),:);
    
    vec_one = matWVLST(tmp_verts(1),:) - tmp_mp;
    vec_two = matWVLST(tmp_verts(2),:) - tmp_mp;
    
    vert_one = ((dot(vec_one, tmp_spandir, 2)).*tmp_spandir) + tmp_mp;
    vert_two = ((dot(vec_two, tmp_spandir, 2)).*tmp_spandir) + tmp_mp;
    
    % Making tmp variables
    tmpWVLST(tmp_verts(1),:) = vert_one;
    tmpWVLST(tmp_verts(2),:) = vert_two;
    
    tmpWCENTER = (tmpWVLST(matWDVE(:,1),:) + tmpWVLST(matWDVE(:,2),:) + tmpWVLST(matWDVE(:,3),:))./3;
    
    % matWPLEX, matWDVECT, matWROTANG
    P = permute(reshape(tmpWVLST(matWDVE(:,:)',:)', 3, 3, []), [2 1 3]);
    DNORM = cross(tmpWVLST(matWDVE(:,2),:) - tmpWVLST(matWDVE(:,3),:), tmpWVLST(matWDVE(:,1),:) - tmpWVLST(matWDVE(:,3),:), 2);
    DNORM = DNORM./sqrt(sum(DNORM.^2,2));
    DNORM(vecWDVEFLIP,:) = DNORM(vecWDVEFLIP,:).*-1;
    [tmpWPLEX, ~, tmpWROTANG] = fcnTRITOLEX(P, DNORM, tmpWCENTER);
    tmpWCOEFF = fcnADJCOEFF(vecWVMU, vecWEMU, tmpWVLST, tmpWCENTER, tmpWROTANG, matWDVE, matWCOEFF, matWELST, matWEIDX, valWNELE);
    
    fpg_og(:,:,i) = locations.*vert_one + (1 - locations).*vert_two;
    wind(:,:,i) = fcnSDVEVEL(fpg_og(:,:,i), valWNELE, tmpWCOEFF, tmpWPLEX, tmpWROTANG, tmpWCENTER, vecWDVESYM, boundind, 1e-8);

end


% fpg2 = reshape(permute(fpg_og, [2 1 3]), size(fpg_og, 2), [])';
% wind2 = reshape(permute(wind, [2 1 3]), size(wind, 2), [])';
% hold on
% quiver3(fpg2(:,1), fpg2(:,2), fpg2(:,3), wind2(:,1), wind2(:,2), wind2(:,3))
% hold off

%
% TAU?????? need to change?
%

% Velocity along the TE
u = []; v = []; w = [];
for i = 1:size(vecTEDVE,1)
    uvw = fcnGLOBSTAR(wind(:,:,i), repmat(matROTANG(vecTEDVE(i),:), length(locations), 1));
    u(i,:) = polyfit(tau(:,1,i), uvw(:,1), 2);
    v(i,:) = polyfit(tau(:,1,i), uvw(:,2), 2);
    w(i,:) = polyfit(tau(:,1,i), uvw(:,3), 2);
end

% Induced lift and side force (Wind x spanwise_direction)
tmp = fcnGLOBSTAR(matSPANDIR(vecTEDVE,:), matROTANG(vecTEDVE,:));
u_out = fcnKJU(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY, tmp(:,1), tmp(:,2), tmp(:,3), lim);
v_out = fcnKJV(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY, tmp(:,1), tmp(:,2), tmp(:,3), lim);
w_out = fcnKJW(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, valDENSITY, tmp(:,1), tmp(:,2), tmp(:,3), lim);
F = fcnSTARGLOB([u_out v_out w_out], matROTANG(vecTEDVE,:));

tmp = matVLST(matELST(vecTE,2),:) - matVLST(matELST(vecTE,1),:);
tmp = abs(sqrt(sum(tmp.^2,2))./(xi_3 - xi_1));

% Splitting special force into lift and side forces
liftind(vecTEDVE,1) = dot(F./tmp, matDVELIFT_DIR(vecTEDVE,:), 2);
sideind(vecTEDVE,1) = dot(F./tmp, matDVESIDE_DIR(vecTEDVE,:), 2);
dragind(vecTEDVE,1) = dot(F./tmp, matDVEDRAG_DIR(vecTEDVE,:), 2);



%% Combining forces
% Apparent mass

if strcmpi(strWAKE_TYPE, 'UNSTEADY')
    
    
    
    
    
end


vecDVELIFT = liftfree + liftind;
vecDVESIDE = sidefree + sideind;
vecDVEDRAG = dragind;

% Symmetry
vecDVELIFT(vecDVESYM) = vecDVELIFT(vecDVESYM)*2;
vecDVEDRAG(vecDVESYM) = vecDVEDRAG(vecDVESYM)*2;

end