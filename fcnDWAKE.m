function [matWCOEFF] = fcnDWAKE(type, valTIMESTEP, strATYPE, vecULS, valWNELE, vecWLE, vecWLEDVE, vecWTE, vecWTEDVE, matWEATT, matWELST, matWROTANG, matWCENTER, matWVLST, vecTE, vecTEDVE, matCOEFF, matCENTER, matROTANG, matWCOEFF, matWPLEX, vecWDVECIRC, vecWSYM, vecWSYMDVE, vecWDVESYM, matWVATT, matWDVE)

 %% Circulation equations between elements
% Evaluated at the mid-point of each edge which splits two HDVEs
idx = all(matWEATT,2); % All edges that split 2 DVEs

vnum_a = matWVLST(matWELST(idx,1),:);
vnum_b = matWVLST(matWELST(idx,2),:);
vnum_mid = (vnum_a + vnum_b)./2;

%% Circulation at edge corner and midpoints
dvenum = [matWEATT(idx,1) matWEATT(idx,2)];
circ = [fcnDCIRC(repmat(vnum_a,1,1,2), dvenum, valWNELE, matWROTANG, matWCENTER); ...
        fcnDCIRC(repmat(vnum_mid,1,1,2), dvenum, valWNELE, matWROTANG, matWCENTER); ...
        fcnDCIRC(repmat(vnum_b,1,1,2), dvenum, valWNELE, matWROTANG, matWCENTER)];
res_circ = zeros(size(circ,1),1);

%% Vorticity along edge between elements
% Unit vector in local ref frame (a for HDVE1, b for HDVE2) from local vertex to local vertex on the edge that forms the border between the two
% vort = [fcnDVORTEDGE(repmat(vnum_a,1,1,2), dvenum, valWNELE, matWROTANG, matWCENTER); ...
%         fcnDVORTEDGE(repmat(vnum_mid,1,1,2), dvenum, valWNELE, matWROTANG, matWCENTER); ...
%         fcnDVORTEDGE(repmat(vnum_b,1,1,2), dvenum, valWNELE, matWROTANG, matWCENTER)];
vort = fcnDVORTEDGE(repmat(vnum_mid,1,1,2), dvenum, valWNELE, matWROTANG, matWCENTER);
res_vort = zeros(size(vort,1),1);
    
%% Leading edge of wake
pts(:,:,1) = matWVLST(matWELST(vecWLE,1),:);
pts(:,:,2) = matWVLST(matWELST(vecWLE,2),:);
pts(:,:,3) = (pts(:,:,1) + pts(:,:,2))./2;

circ_le = [   fcnDCIRC2(pts(:,:,1), vecWLEDVE, valWNELE, matWROTANG, matWCENTER); ...
    fcnDCIRC2(pts(:,:,2), vecWLEDVE, valWNELE, matWROTANG, matWCENTER); ...
    fcnDCIRC2(pts(:,:,3), vecWLEDVE, valWNELE, matWROTANG, matWCENTER)];

pts_loc(:,:,1) = fcnGLOBSTAR(pts(:,:,1) - matCENTER(vecTEDVE(:,1),:), matROTANG(vecTEDVE(:,1),:));
pts_loc(:,:,2) = fcnGLOBSTAR(pts(:,:,2) - matCENTER(vecTEDVE(:,1),:), matROTANG(vecTEDVE(:,1),:));
pts_loc(:,:,3) = fcnGLOBSTAR(pts(:,:,3) - matCENTER(vecTEDVE(:,1),:), matROTANG(vecTEDVE(:,1),:));
res1_1 = [   sum([0.5.*pts_loc(:,2,1).^2 pts_loc(:,2,1) 0.5.*pts_loc(:,1,1).^2 pts_loc(:,1,1) pts_loc(:,2,1).*pts_loc(:,1,1) ones(size(pts_loc(:,1,1)))].*matCOEFF(vecTEDVE(:,1),:),2); ...
    sum([0.5.*pts_loc(:,2,2).^2 pts_loc(:,2,2) 0.5.*pts_loc(:,1,2).^2 pts_loc(:,1,2) pts_loc(:,2,2).*pts_loc(:,1,2) ones(size(pts_loc(:,1,2)))].*matCOEFF(vecTEDVE(:,1),:),2); ...
    sum([0.5.*pts_loc(:,2,3).^2 pts_loc(:,2,3) 0.5.*pts_loc(:,1,3).^2 pts_loc(:,1,3) pts_loc(:,2,3).*pts_loc(:,1,3) ones(size(pts_loc(:,1,3)))].*matCOEFF(vecTEDVE(:,1),:),2)];
if strcmpi(strATYPE{2},'PANEL')
    pts_loc(:,:,1) = fcnGLOBSTAR(pts(:,:,1) - matCENTER(vecTEDVE(:,2),:), matROTANG(vecTEDVE(:,2),:));
    pts_loc(:,:,2) = fcnGLOBSTAR(pts(:,:,2) - matCENTER(vecTEDVE(:,2),:), matROTANG(vecTEDVE(:,2),:));
    pts_loc(:,:,3) = fcnGLOBSTAR(pts(:,:,3) - matCENTER(vecTEDVE(:,2),:), matROTANG(vecTEDVE(:,2),:));
    res1_2 = [   sum([0.5.*pts_loc(:,2,1).^2 pts_loc(:,2,1) 0.5.*pts_loc(:,1,1).^2 pts_loc(:,1,1) pts_loc(:,2,1).*pts_loc(:,1,1) ones(size(pts_loc(:,1,1)))].*matCOEFF(vecTEDVE(:,2),:),2); ...
        sum([0.5.*pts_loc(:,2,2).^2 pts_loc(:,2,2) 0.5.*pts_loc(:,1,2).^2 pts_loc(:,1,2) pts_loc(:,2,2).*pts_loc(:,1,2) ones(size(pts_loc(:,1,2)))].*matCOEFF(vecTEDVE(:,2),:),2); ...
        sum([0.5.*pts_loc(:,2,3).^2 pts_loc(:,2,3) 0.5.*pts_loc(:,1,3).^2 pts_loc(:,1,3) pts_loc(:,2,3).*pts_loc(:,1,3) ones(size(pts_loc(:,1,3)))].*matCOEFF(vecTEDVE(:,2),:),2)];
    
    res_circ_le = res1_2 + res1_1;
else
    res_circ_le = res1_1;
end

% vort_le = [fcnDVORT2(pts(:,:,1), vecWLEDVE, valWNELE, matWCENTER, matWROTANG, 'A');...
%     fcnDVORT2(pts(:,:,2), vecWLEDVE, valWNELE, matWCENTER, matWROTANG, 'A');...
%     fcnDVORT2(pts(:,:,3), vecWLEDVE, valWNELE, matWCENTER, matWROTANG, 'A')];
vort_le = fcnDVORT2(pts(:,:,2), vecWLEDVE, valWNELE, matWCENTER, matWROTANG, 'A');
res_vort_le = zeros(size(vort_le,1),1);
  

%% Integrated circulation
circ_int = [];
res_circ_int = [];
vort_steady = [];
res_steady = [];
if strcmpi(type, 'UNSTEADY')
    xi_1 = permute(matWPLEX(1,1,:),[3 2 1]);
    xi_2 = permute(matWPLEX(2,1,:),[3 2 1]);
    xi_3 = permute(matWPLEX(3,1,:),[3 2 1]);

    eta_1 = permute(matWPLEX(1,2,:),[3 2 1]);
    eta_2 = permute(matWPLEX(2,2,:),[3 2 1]);
    eta_3 = permute(matWPLEX(3,2,:),[3 2 1]);

    idx_flp = abs(xi_2 - xi_3) < 1e-5;
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

    a1 = -((xi_1-xi_3).*((C.^3-E.^3).*xi_1.^3+((C.^3-E.^3).*xi_3+4.*C.^2.*D_LE-4.*E.^2.*D_TE).*xi_1.^2+((C.^3-E.^3).*xi_3.^2+(4.*C.^2.*D_LE-4.*E.^2.*D_TE).*xi_3+6.*C.*D_LE.^2-6.*E.*D_TE.^2).*xi_1+(C.^3-E.^3).*xi_3.^3+(4.*C.^2.*D_LE-4.*E.^2.*D_TE).*xi_3.^2+(6.*C.*D_LE.^2-6.*E.*D_TE.^2).*xi_3+4.*D_LE.^3-4.*D_TE.^3))./0.24e2;
    a2 = -((xi_1-xi_3).*((C.^2-E.^2).*xi_1.^2+((C.^2-E.^2).*xi_3+3.*C.*D_LE-3.*E.*D_TE).*xi_1+(C.^2-E.^2).*xi_3.^2+(3.*C.*D_LE-3.*E.*D_TE).*xi_3+3.*D_LE.^2-3.*D_TE.^2))./0.6e1;
    b1 = -(xi_1-xi_3).*((C-E).*xi_1.^3+((C-E).*xi_3+0.4e1./0.3e1.*D_LE-0.4e1./0.3e1.*D_TE).*xi_1.^2+((C-E).*xi_3+0.4e1./0.3e1.*D_LE-0.4e1./0.3e1.*D_TE).*xi_3.*xi_1+((C-E).*xi_3+0.4e1./0.3e1.*D_LE-0.4e1./0.3e1.*D_TE).*xi_3.^2)./0.8e1;
    b2 = -(xi_1-xi_3).*((C-E).*xi_1.^2+((C-E).*xi_3+0.3e1./0.2e1.*D_LE-0.3e1./0.2e1.*D_TE).*xi_1+((C-E).*xi_3+0.3e1./0.2e1.*D_LE-0.3e1./0.2e1.*D_TE).*xi_3)./0.3e1;
    c2 = -(xi_1-xi_3).*((C.^2-E.^2).*xi_1.^3+((C.^2-E.^2).*xi_3+0.8e1./0.3e1.*C.*D_LE-0.8e1./0.3e1.*E.*D_TE).*xi_1.^2+((C.^2-E.^2).*xi_3.^2+(0.8e1./0.3e1.*C.*D_LE-0.8e1./0.3e1.*E.*D_TE).*xi_3+0.2e1.*D_LE.^2-0.2e1.*D_TE.^2).*xi_1+xi_3.*((C.^2-E.^2).*xi_3.^2+(0.8e1./0.3e1.*C.*D_LE-0.8e1./0.3e1.*E.*D_TE).*xi_3+0.2e1.*D_LE.^2-0.2e1.*D_TE.^2))./0.8e1;
    c3 = -((xi_1-xi_3).*((C-E).*xi_1+(C-E).*xi_3+2.*D_LE-2.*D_TE))./0.2e1;
    gamma = [a1 a2 b1 b2 c2 c3];

    circ_int = fcnCREATEDSECT(sparse(size(gamma,1), valWNELE*6), size(gamma,1), 6, [1:valWNELE]', [], gamma, []);
    res_circ_int = vecWDVECIRC;
    res_circ_int(idx_flp,:) = res_circ_int(idx_flp,:).*-1; 
else
    vort_steady = fcnDVORT1([1:valWNELE]', valWNELE, 'A');
    res_steady = zeros(size(vort_steady,1),1); 
end

%%
DW = [circ; vort; circ_le; vort_le; circ_int; vort_steady];
RW = [res_circ; res_vort; res_circ_le; res_vort_le; res_circ_int; res_steady];

matWCOEFF = DW\RW;
matWCOEFF = reshape(matWCOEFF,6,valWNELE,1)';

end

