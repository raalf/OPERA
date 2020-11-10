function [v_xa, v_xb, v_xc, v_ea, v_eb, v_ec, v_za, v_zb, v_zc] = fcnVMAP(valNELE, u, matKINCON_DVE, matKINCON_P, matCENTER, matROTANG)
%For approximating the surface velocity distribution across a DDE using
% a first order approx in both the spanwise and chordwise directions

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