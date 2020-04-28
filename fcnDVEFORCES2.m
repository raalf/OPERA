function [matF_FS, matF_IF, matF_ID, matINTCIRC] = fcnDVEFORCES2(valTIMESTEP, matVLST, matELST, matROTANG, ...
    matCOEFF, vecTEDVE, valDENSITY, valNELE, matSPANDIR, vecTE, matPLEX, matWCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG, ...
    matWVLST, vecWLEDVE, matWELST, vecDVESYM, vecWDVESYM, valWSIZE, matWDVE, vecWDVEFLIP, matWEIDX, vecWEMU, ...
    vecWVMU, matCENTER, matKINCON_P, matKINCON_DVE, matUINF_KK, vecWLE, matWDVECT)

%% Common params
xi_1 = permute(matPLEX(1,1,:),[3 2 1]);
xi_2 = permute(matPLEX(2,1,:),[3 2 1]);
xi_3 = permute(matPLEX(3,1,:),[3 2 1]);

eta_1 = permute(matPLEX(1,2,:),[3 2 1]);
eta_2 = permute(matPLEX(2,2,:),[3 2 1]);
eta_3 = permute(matPLEX(3,2,:),[3 2 1]);

idx_flp = xi_3 < xi_1; % Flipping influence of elements that need a good flippin

C = (eta_3 - eta_2)./(xi_3 - xi_2);
D_LE = eta_2 - ((xi_2.*(eta_3 - eta_2))./(xi_3 - xi_2));
E = (eta_3 - eta_1)./(xi_3 - xi_1);
D_TE = eta_1 - ((xi_1.*(eta_3 - eta_1))./(xi_3 - xi_1));

A_1 = matCOEFF(:,1); A_2 = matCOEFF(:,2); B_1 = matCOEFF(:,3);
B_2 = matCOEFF(:,4); C_2 = matCOEFF(:,5); C_3 = matCOEFF(:,6);

len_adj = matVLST(matELST(vecTE,2),:) - matVLST(matELST(vecTE,1),:);
len_adj = abs(sqrt(sum(len_adj.^2,2))./(xi_3(vecTEDVE) - xi_1(vecTEDVE)));

%% Freestream and Induced Lift and Sideforce
% Mapping freestream and surface+wake induced velocities across element
% surfaces. 9 unknowns for velocity across each element (3 per direction)
% U_ES = (eta*v__xb + v__xa*xi + v__xc, eta*v__eb + v__ea*xi + v__ec, eta*v__zb + v__za*xi + v__zc)

matF_FS = fcnDVEFORCE(idx_flp, valNELE, valDENSITY, matUINF_KK, matKINCON_DVE, matKINCON_P, matCENTER, matROTANG, A_1, A_2, B_1, B_2, C_2, xi_1, xi_3, C, D_LE, E, D_TE);

tmp_w = fcnSDVEVEL(matKINCON_P, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM, [], 0) + fcnSDVEVEL(matKINCON_P, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, vecWDVESYM, [], 0);
matF_IF = fcnDVEFORCE(idx_flp, valNELE, valDENSITY, tmp_w, matKINCON_DVE, matKINCON_P, matCENTER, matROTANG, A_1, A_2, B_1, B_2, C_2, xi_1, xi_3, C, D_LE, E, D_TE);

%% Induced Drag
% Induced velocities at wake leading edge DVEs (wind is 3 x 3 x num_dve)

locations = linspace(0.2, 0.8, 3)';

wind = nan(size(locations,1), 3, size(vecWLEDVE,1));
fpg_og = nan(size(locations,1), 3, size(vecWLEDVE,1));
fpg_us = nan(size(locations,1), 3, size(vecWLEDVE,1));

boundind = false(valWNELE,1);
boundind(vecWLEDVE) = true;

% Need to adjust wake LE elements to straight, one by one
offset = 1e-2;
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
    [tmpWPLEX, tmpWDVECT, tmpWROTANG] = fcnTRITOLEX(P, DNORM, tmpWCENTER);
    tmpWCOEFF = fcnADJCOEFF(vecWVMU, vecWEMU, tmpWVLST, tmpWCENTER, tmpWROTANG, matWDVE, matWCOEFF, matWELST, matWEIDX, valWNELE);
    
    fpg_og(:,:,i) = locations.*matWVLST(tmp_verts(1),:) + (1 - locations).*(matWVLST(tmp_verts(2),:));
    fpg_us(:,:,i) = locations.*vert_one + (1 - locations).*vert_two;
    wind(:,:,i) = (fcnSDVEVEL(fpg_us(:,:,i) + tmpWDVECT(vecWLEDVE(i),:,3).*offset, valWNELE, tmpWCOEFF, tmpWPLEX, tmpWROTANG, tmpWCENTER, vecWDVESYM, boundind, 0) + ...
        fcnSDVEVEL(fpg_us(:,:,i) - tmpWDVECT(vecWLEDVE(i),:,3).*offset, valWNELE, tmpWCOEFF, tmpWPLEX, tmpWROTANG, tmpWCENTER, vecWDVESYM, boundind, 0))./2;
end

% fpg2 = reshape(permute(fpg_og, [2 1 3]), size(fpg_og, 2), [])';
% wind2 = reshape(permute(wind, [2 1 3]), size(wind, 2), [])';
% hold on
% quiver3(fpg2(:,1), fpg2(:,2), fpg2(:,3), wind2(:,1), wind2(:,2), wind2(:,3))
% hold off
% view([90 0])

% Velocity along the TE
u = []; v = []; w = [];
for i = 1:size(vecTEDVE,1)
    fpl = fcnGLOBSTAR(fpg_og(:,:,i) - matCENTER(vecTEDVE(i),:), matROTANG(vecTEDVE(i),:));
    uvw = fcnGLOBSTAR(wind(:,:,i), repmat(matROTANG(vecTEDVE(i),:), length(locations), 1));
    u(i,:) = polyfit(fpl(:,1), uvw(:,1), 2);
    v(i,:) = polyfit(fpl(:,1), uvw(:,2), 2);
    w(i,:) = polyfit(fpl(:,1), uvw(:,3), 2);
end

% Induced lift and side force (Wind x spanwise_direction)
lim = 1e10;
tmp = fcnGLOBSTAR(matSPANDIR(vecTEDVE,:), matROTANG(vecTEDVE,:));
u_out = fcnKJU(xi_1(vecTEDVE), xi_3(vecTEDVE), eta_1(vecTEDVE), eta_3(vecTEDVE), E(vecTEDVE), D_TE(vecTEDVE), A_1(vecTEDVE), A_2(vecTEDVE), B_1(vecTEDVE), B_2(vecTEDVE), C_2(vecTEDVE), C_3(vecTEDVE), u, v, w, valDENSITY, tmp(:,1), tmp(:,2), tmp(:,3), lim);
v_out = fcnKJV(xi_1(vecTEDVE), xi_3(vecTEDVE), eta_1(vecTEDVE), eta_3(vecTEDVE), E(vecTEDVE), D_TE(vecTEDVE), A_1(vecTEDVE), A_2(vecTEDVE), B_1(vecTEDVE), B_2(vecTEDVE), C_2(vecTEDVE), C_3(vecTEDVE), u, v, w, valDENSITY, tmp(:,1), tmp(:,2), tmp(:,3), lim);
w_out = fcnKJW(xi_1(vecTEDVE), xi_3(vecTEDVE), eta_1(vecTEDVE), eta_3(vecTEDVE), E(vecTEDVE), D_TE(vecTEDVE), A_1(vecTEDVE), A_2(vecTEDVE), B_1(vecTEDVE), B_2(vecTEDVE), C_2(vecTEDVE), C_3(vecTEDVE), u, v, w, valDENSITY, tmp(:,1), tmp(:,2), tmp(:,3), lim);
F = fcnSTARGLOB([u_out v_out w_out], matROTANG(vecTEDVE,:));

% Splitting special force into lift and side forces
matF_ID(vecTEDVE,:) = F./len_adj;

%% Integrated circulation
tmp = fcnINTCIRC(xi_1, xi_3, C, D_LE, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3);
tmp(~idx_flp) = -tmp(~idx_flp);
matINTCIRC = tmp;

end

function [v_xa, v_xb, v_xc, v_ea, v_eb, v_ec, v_za, v_zb, v_zc] = fcnVMAP(valNELE, u, matKINCON_DVE, matKINCON_P, matCENTER, matROTANG)
    % points and their velocities
    kk_loc = fcnGLOBSTAR(matKINCON_P - matCENTER(matKINCON_DVE,:), matROTANG(matKINCON_DVE,:));
    uinf_loc = fcnGLOBSTAR(u, matROTANG(matKINCON_DVE,:));

    v_xa = nan(valNELE,1); v_xb = nan(valNELE,1); v_xc = nan(valNELE,1);
    v_ea = nan(valNELE,1); v_eb = nan(valNELE,1); v_ec = nan(valNELE,1);
    v_za = nan(valNELE,1); v_zb = nan(valNELE,1); v_zc = nan(valNELE,1);
    for i = 1:valNELE
        idx = matKINCON_DVE == i;

        pts = kk_loc(idx,:);
        pts(:,3) = 1;
        vs = uinf_loc(idx,:);

        vmap_fs = pts \ vs;

        v_xa(i,1) = vmap_fs(1,1);
        v_xb(i,1) = vmap_fs(2,1);
        v_xc(i,1) = vmap_fs(3,1);
        v_ea(i,1) = vmap_fs(1,2);
        v_eb(i,1) = vmap_fs(2,2);
        v_ec(i,1) = vmap_fs(3,2);
        v_za(i,1) = vmap_fs(1,3);
        v_zb(i,1) = vmap_fs(2,3);
        v_zc(i,1) = vmap_fs(3,3);
    end
end

function f_out = fcnDVEFORCE(idx_flp, valNELE, valDENSITY, u, matKINCON_DVE, matKINCON_P, matCENTER, matROTANG, A_1, A_2, B_1, B_2, C_2, xi_1, xi_3, C, D_LE, E, D_TE)
    [v_xa, v_xb, v_xc, v_ea, v_eb, v_ec, v_za, v_zb, v_zc] = fcnVMAP(valNELE, u, matKINCON_DVE, matKINCON_P, matCENTER, matROTANG);

    t1 = (valDENSITY .* v_zb);
    t2 = (C .^ 2);
    t4 = (E .^ 2);
    t10 = valDENSITY .* B_1;
    t12 = valDENSITY .* C_2;
    t14 = v_zb .* t10 + v_za .* t12;
    t15 = -t2 + t4;
    t18 = valDENSITY .* v_za;
    t19 = -C + E;
    t23 = xi_1 .^ 2;
    t24 = t23 .^ 2;
    t25 = xi_3 .^ 2;
    t26 = t25 .^ 2;
    t37 = valDENSITY .* B_2;
    t40 = v_zc .* t12 + v_zb .* t37;
    t46 = -2 .* C .* D_LE + 2 .* D_TE .* E;
    t51 = -v_zc .* t10 - v_za .* t37;
    t53 = -D_LE + D_TE;
    t62 = D_LE .^ 2;
    t64 = D_TE .^ 2;
    t73 = -t62 + t64;
    t76 = valDENSITY .* v_zc;
    t88 = (xi_3 - xi_1);
    fx = (-t24 + t26) .* (((-t2 .* C + t4 .* E) .* C_2 .* t1) ./ 0.3e1 + t15 .* t14 ./ 0.2e1 + t19 .* B_1 .* t18) ./ ...
        0.4e1 + (-t23 .* xi_1 + t25 .* xi_3) .* (((-3 .* D_LE .* t2 + 3 .* t4 .* D_TE) .* C_2 .* t1) ./ 0.3e1 + t15 .* ...
        t40 ./ 0.2e1 + t46 .* t14 ./ 0.2e1 - t19 .* t51 + t53 .* B_1 .* t18) ./ 0.3e1 + (-t23 + t25) .* (((-3 .* t62 .* C + ...
        3 .* E .* t64) .* C_2 .* t1) ./ 0.3e1 + t46 .* t40 ./ 0.2e1 + t73 .* t14 ./ 0.2e1 + t19 .* B_2 .* t76 - t53 .* t51) ...
        ./ 0.2e1 + (t88 .* (-t62 .* D_LE + t64 .* D_TE) .* C_2 .* t1) ./ 0.3e1 + t88 .* t73 .* t40 ./ 0.2e1 + t88 .* t53 .* B_2 .* t76;

    t1 = (valDENSITY .* v_zb);
    t2 = (C .^ 2);
    t4 = (E .^ 2);
    t10 = valDENSITY .* A_1;
    t12 = valDENSITY .* C_2;
    t14 = v_za .* t10 + v_zb .* t12;
    t15 = -t2 + t4;
    t18 = valDENSITY .* v_za;
    t19 = -C + E;
    t23 = xi_1 .^ 2;
    t24 = t23 .^ 2;
    t25 = xi_3 .^ 2;
    t26 = t25 .^ 2;
    t38 = valDENSITY .* A_2;
    t40 = v_zc .* t10 + v_zb .* t38;
    t46 = -2 .* C .* D_LE + 2 .* D_TE .* E;
    t51 = v_zc .* t12 + v_za .* t38;
    t53 = -D_LE + D_TE;
    t62 = D_LE .^ 2;
    t64 = D_TE .^ 2;
    t73 = -t62 + t64;
    t76 = valDENSITY .* v_zc;
    t88 = (xi_3 - xi_1);
    fe = (-t24 + t26) .* (((-t2 .* C + t4 .* E) .* A_1 .* t1) ./ 0.3e1 + t15 .* t14 ./ 0.2e1 + t19 .* C_2 .* t18) ./ 0.4e1 + ...
        (-t23 .* xi_1 + t25 .* xi_3) .* (((-3 .* D_LE .* t2 + 3 .* t4 .* D_TE) .* A_1 .* t1) ./ 0.3e1 + t15 .* t40 ./ 0.2e1 + t46 .* ...
        t14 ./ 0.2e1 + t19 .* t51 + t53 .* C_2 .* t18) ./ 0.3e1 + (-t23 + t25) .* (((-3 .* t62 .* C + 3 .* E .* t64) .* A_1 .* t1) ./ ...
        0.3e1 + t46 .* t40 ./ 0.2e1 + t73 .* t14 ./ 0.2e1 + t19 .* A_2 .* t76 + t53 .* t51) ./ 0.2e1 + (t88 .* (-t62 .* D_LE + t64 .*...
        D_TE) .* A_1 .* t1) ./ 0.3e1 + t88 .* t73 .* t40 ./ 0.2e1 + t88 .* t53 .* A_2 .* t76;

    t1 = valDENSITY .* v_eb;
    t3 = valDENSITY .* v_xb;
    t5 = (-A_1 .* t1 - C_2 .* t3);
    t6 = (C .^ 2);
    t8 = (E .^ 2);
    t13 = valDENSITY .* v_ea;
    t17 = valDENSITY .* v_xa;
    t19 = -A_1 .* t13 - B_1 .* t3 - C_2 .* t1 - C_2 .* t17;
    t20 = -t6 + t8;
    t23 = -C + E;
    t29 = xi_1 .^ 2;
    t30 = t29 .^ 2;
    t31 = xi_3 .^ 2;
    t32 = t31 .^ 2;
    t42 = valDENSITY .* v_ec;
    t46 = valDENSITY .* v_xc;
    t48 = -A_1 .* t42 - A_2 .* t1 - B_2 .* t3 - C_2 .* t46;
    t54 = -2 .* D_LE .* C + 2 .* D_TE .* E;
    t59 = -B_1 .* t46 - B_2 .* t17;
    t61 = -D_LE + D_TE;
    t66 = A_2 .* t13 + C_2 .* t42;
    t76 = D_LE .^ 2;
    t78 = D_TE .^ 2;
    t86 = -t76 + t78;
    t103 = (xi_3 - xi_1);
    fz = (-t30 + t32) .* (((-t6 .* C + t8 .* E) .* t5) ./ 0.3e1 + t20 .* t19 ./ 0.2e1 - t23 .* B_1 .* t17 - t23 .* C_2 .* t13) ./ 0.4e1 + (-t29 .* xi_1 + t31 .* xi_3) .* ...
        (((-3 .* t6 .* D_LE + 3 .* t8 .* D_TE) .* t5) ./ 0.3e1 + t20 .* t48 ./ 0.2e1 + t54 .* t19 ./ 0.2e1 + t23 .* t59 - t61 .* B_1 .* t17 - t23 .* t66 - t61 .* C_2 .* t13) ./ ...
        0.3e1 + (-t29 + t31) .* (((-3 .* C .* t76 + 3 .* E .* t78) .* t5) ./ 0.3e1 + t54 .* t48 ./ 0.2e1 + t86 .* t19 ./ 0.2e1 - t23 .* B_2 .* t46 + t61 .* t59 - t23 .* A_2 .* ...
        t42 - t61 .* t66) ./ 0.2e1 + (t103 .* (-t76 .* D_LE + t78 .* D_TE) .* t5) ./ 0.3e1 + t103 .* t86 .* t48 ./ 0.2e1 - t103 .* t61 .* B_2 .* t46 - t103 .* t61 .* A_2 .* t42;

    f_out = [fx fe fz];
f_out(idx_flp,:) = f_out(idx_flp,:).*-1;

f_out = fcnSTARGLOB(f_out, matROTANG);

end

















