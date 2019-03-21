function [matWCOEFF, vecWDVECIRCnew] = fcnDWAKENEW(valTIMESTEP, strATYPE, vecULS, valWNELE, vecWLE, vecWLEDVE, vecWTE, vecWTEDVE, matWEATT, matWELST, matWROTANG, matWCENTER, matWVLST, vecTE, vecTEDVE, matCOEFF, matCENTER, matROTANG, matWCOEFF, matWPLEX, vecWSYM, vecWSYMDVE, vecWDVESYM)

valWNELE_tmp = length([vecWLEDVE; vecWTEDVE]);

%% Spanwise circulation in the front half of newest wake row
% Setting circulation at LE of wake equal to TE of wing

% tmpvrt = fcnDVORT1(vecWLEDVE, valWNELE, 'B');
% nedg = size(vecWLEDVE,1);
% zer = zeros(nedg,1);
% c_term = [zer, zer, zer, zer, ones(nedg,1)];
% tmpvrt2 = fcnCREATEDSECT(sparse(nedg, valWNELE*5), nedg, 5, vecWLEDVE, [], c_term, []);
% vort1 = [tmpvrt; tmpvrt2];
% res1 = -[matCOEFF(vecTEDVE,3); matCOEFF(vecTEDVE,4); matCOEFF(vecTEDVE,5)];

pts(:,:,1) = matWVLST(matWELST(vecWLE,1),:);
pts(:,:,2) = matWVLST(matWELST(vecWLE,2),:);
pts(:,:,3) = (pts(:,:,1) + pts(:,:,2))./2;

circ1 = [   fcnDCIRC2(pts(:,:,1), vecWLEDVE, valWNELE, matWROTANG, matWCENTER); ...
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
    
    res1 = res1_2 + res1_1;
else
    res1 = res1_1;
end

vort6 = [fcnDVORT2(pts(:,:,1), vecWLEDVE, valWNELE, matWCENTER, matWROTANG, 'A');...
    fcnDVORT2(pts(:,:,2), vecWLEDVE, valWNELE, matWCENTER, matWROTANG, 'A')];
% vort6 = [fcnDVORT2(pts(:,:,1), vecWLEDVE, valWNELE, matWCENTER, matWROTANG, 'A');...
%     fcnDVORT2(pts(:,:,2), vecWLEDVE, valWNELE, matWCENTER, matWROTANG, 'A');...
%     fcnDVORT2(pts(:,:,3), vecWLEDVE, valWNELE, matWCENTER, matWROTANG, 'A')];
res6 = zeros(size(vort6,1),1);

%% Spanwise circulation in the back half of the newest wake row
% Circulation equations between elements
% Evaluated at the mid-point of each edge which splits two HDVEs in most
% recent row of wake DVEs
[~,idx] = ismember(sort([vecWLEDVE vecWTEDVE],2), sort(matWEATT,2),'rows');
vnum_a = matWVLST(matWELST(idx,1),:);
vnum_b = matWVLST(matWELST(idx,2),:);
vnum_mid = (vnum_a + vnum_b)./2;

% Circulation at edge corner and midpoints
dvenum = [matWEATT(idx,1) matWEATT(idx,2)];
circ2 = [   fcnDCIRC(repmat(vnum_a,1,1,2), dvenum, valWNELE, matWROTANG, matWCENTER); ...
    fcnDCIRC(repmat(vnum_mid,1,1,2), dvenum, valWNELE, matWROTANG, matWCENTER); ...
    fcnDCIRC(repmat(vnum_b,1,1,2), dvenum, valWNELE, matWROTANG, matWCENTER)];
res2 = zeros(size(circ2,1),1);

% Vorticity along edge between elements
% Unit vector in local ref frame (a for HDVE1, b for HDVE2) from local vertex to local vertex on the edge that forms the border between the two
vort3 = fcnDVORTEDGE(repmat(vnum_mid,1,1,2), dvenum, valWNELE, matWROTANG, matWCENTER);
% vort3 = [   fcnDVORTEDGE(repmat(vnum_a,1,1,2), dvenum, valWNELE, matWROTANG, matWCENTER); ...
%     fcnDVORTEDGE(repmat(vnum_mid,1,1,2), dvenum, valWNELE, matWROTANG, matWCENTER); ...
%     fcnDVORTEDGE(repmat(vnum_b,1,1,2), dvenum, valWNELE, matWROTANG, matWCENTER)];
res3 = zeros(size(vort3,1),1);


%% Chordwise circulation in the newest wake row
vort4 = [];
res4 = [];
if valTIMESTEP == 1
    dvenum = [vecWLEDVE; vecWTEDVE];
    vort4 = fcnDVORT1(dvenum, valWNELE, 'A');
    res4 = zeros(size(vort4,1),1);
    
    circ5 = [];
    res5 = [];
else
    % Set circulation at TE of newest wake row equal to LE of 2nd newest
    % wake row
    
    pts(:,:,1) = matWVLST(matWELST(vecWTE,1),:);
    pts(:,:,2) = matWVLST(matWELST(vecWTE,2),:);
    pts(:,:,3) = (pts(:,:,1) + pts(:,:,2))./2;
    
    circ5 = [   fcnDCIRC2(pts(:,:,1), vecWTEDVE, valWNELE, matWROTANG, matWCENTER); ...
        fcnDCIRC2(pts(:,:,2), vecWTEDVE, valWNELE, matWROTANG, matWCENTER); ...
        fcnDCIRC2(pts(:,:,3), vecWTEDVE, valWNELE, matWROTANG, matWCENTER)];
    
    dvenum = vecWTEDVE - length(vecWTEDVE)*3;
    pts_loc(:,:,1) = fcnGLOBSTAR(pts(:,:,1) - matWCENTER(dvenum,:), matWROTANG(dvenum,:));
    pts_loc(:,:,2) = fcnGLOBSTAR(pts(:,:,2) - matWCENTER(dvenum,:), matWROTANG(dvenum,:));
    pts_loc(:,:,3) = fcnGLOBSTAR(pts(:,:,3) - matWCENTER(dvenum,:), matWROTANG(dvenum,:));
    res5 =  [   sum([0.5.*pts_loc(:,2,1).^2 pts_loc(:,2,1) 0.5.*pts_loc(:,1,1).^2 pts_loc(:,1,1) pts_loc(:,2,1).*pts_loc(:,1,1) ones(size(pts_loc(:,1,1)))].*matWCOEFF(dvenum,:),2); ...
        sum([0.5.*pts_loc(:,2,2).^2 pts_loc(:,2,2) 0.5.*pts_loc(:,1,2).^2 pts_loc(:,1,2) pts_loc(:,2,2).*pts_loc(:,1,2) ones(size(pts_loc(:,1,2)))].*matWCOEFF(dvenum,:),2); ...
        sum([0.5.*pts_loc(:,2,3).^2 pts_loc(:,2,3) 0.5.*pts_loc(:,1,3).^2 pts_loc(:,1,3) pts_loc(:,2,3).*pts_loc(:,1,3) ones(size(pts_loc(:,1,3)))].*matWCOEFF(dvenum,:),2)];
    
    %     vort4 = [fcnDVORT2(pts(:,:,1), vecWTEDVE, valWNELE, matWCENTER, matWROTANG, 'A');...
    %         fcnDVORT2(pts(:,:,2), vecWTEDVE, valWNELE, matWCENTER, matWROTANG, 'A');...
    %         fcnDVORT2(pts(:,:,3), vecWTEDVE, valWNELE, matWCENTER, matWROTANG, 'A')];
    %     res4 = [sum([pts_loc(:,2,1) ones(size(pts_loc(:,2,1),1),1), zeros(size(pts_loc(:,2,1),1),3)].*matWCOEFF(dvenum,:),2);...
    %         sum([pts_loc(:,2,2) ones(size(pts_loc(:,2,2),1),1), zeros(size(pts_loc(:,2,2),1),3)].*matWCOEFF(dvenum,:),2);...
    %         sum([pts_loc(:,2,3) ones(size(pts_loc(:,2,3),1),1), zeros(size(pts_loc(:,2,3),1),3)].*matWCOEFF(dvenum,:),2)];
    
    %     nedg = size([vecWLEDVE; vecWTEDVE],1);
    %     zer = zeros(nedg,1);
    %     vort4 = fcnCREATEDSECT(sparse(nedg, valWNELE*5), nedg, 5, [vecWLEDVE; vecWTEDVE], [], [ones(nedg,1), zer, zer, zer, zer], []);
    %     res4 = zeros(size(vort4,1),1);
    
    
end

%% Solving
DW = [circ1; circ2; vort3; vort4; circ5; vort6];
RW = [res1; res2; res3; res4; res5; res6];

DW = DW(:,((min([vecWLEDVE; vecWTEDVE]) - 1)*6 + 1):end);

matWCOEFF = DW\RW;
matWCOEFF = reshape(matWCOEFF,6,valWNELE_tmp,1)';

%% Wake element integrated strength
dves = [vecWLEDVE; vecWTEDVE];

xi_1 = permute(matWPLEX(1,1,dves),[3 2 1]);
xi_2 = permute(matWPLEX(2,1,dves),[3 2 1]);
xi_3 = permute(matWPLEX(3,1,dves),[3 2 1]);

eta_1 = permute(matWPLEX(1,2,dves),[3 2 1]);
eta_2 = permute(matWPLEX(2,2,dves),[3 2 1]);
eta_3 = permute(matWPLEX(3,2,dves),[3 2 1]);

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

a1 = -((xi_1-xi_3).*((C.^3-E.^3).*xi_1.^3+((C.^3-E.^3).*xi_3+4.*C.^2.*D_LE-4.*E.^2.*D_TE).*xi_1.^2+((C.^3-E.^3).*xi_3.^2+(4.*C.^2.*D_LE-4.*E.^2.*D_TE).*xi_3+6.*C.*D_LE.^2-6.*E.*D_TE.^2).*xi_1+(C.^3-E.^3).*xi_3.^3+(4.*C.^2.*D_LE-4.*E.^2.*D_TE).*xi_3.^2+(6.*C.*D_LE.^2-6.*E.*D_TE.^2).*xi_3+4.*D_LE.^3-4.*D_TE.^3))./0.24e2;
a2 = -((xi_1-xi_3).*((C.^2-E.^2).*xi_1.^2+((C.^2-E.^2).*xi_3+3.*C.*D_LE-3.*E.*D_TE).*xi_1+(C.^2-E.^2).*xi_3.^2+(3.*C.*D_LE-3.*E.*D_TE).*xi_3+3.*D_LE.^2-3.*D_TE.^2))./0.6e1;
b1 = -(xi_1-xi_3).*((C-E).*xi_1.^3+((C-E).*xi_3+0.4e1./0.3e1.*D_LE-0.4e1./0.3e1.*D_TE).*xi_1.^2+((C-E).*xi_3+0.4e1./0.3e1.*D_LE-0.4e1./0.3e1.*D_TE).*xi_3.*xi_1+((C-E).*xi_3+0.4e1./0.3e1.*D_LE-0.4e1./0.3e1.*D_TE).*xi_3.^2)./0.8e1;
b2 = -(xi_1-xi_3).*((C-E).*xi_1.^2+((C-E).*xi_3+0.3e1./0.2e1.*D_LE-0.3e1./0.2e1.*D_TE).*xi_1+((C-E).*xi_3+0.3e1./0.2e1.*D_LE-0.3e1./0.2e1.*D_TE).*xi_3)./0.3e1;
c2 = -(xi_1-xi_3).*((C.^2-E.^2).*xi_1.^3+((C.^2-E.^2).*xi_3+0.8e1./0.3e1.*C.*D_LE-0.8e1./0.3e1.*E.*D_TE).*xi_1.^2+((C.^2-E.^2).*xi_3.^2+(0.8e1./0.3e1.*C.*D_LE-0.8e1./0.3e1.*E.*D_TE).*xi_3+0.2e1.*D_LE.^2-0.2e1.*D_TE.^2).*xi_1+xi_3.*((C.^2-E.^2).*xi_3.^2+(0.8e1./0.3e1.*C.*D_LE-0.8e1./0.3e1.*E.*D_TE).*xi_3+0.2e1.*D_LE.^2-0.2e1.*D_TE.^2))./0.8e1;
c3 = -((xi_1-xi_3).*((C-E).*xi_1+(C-E).*xi_3+2.*D_LE-2.*D_TE))./0.2e1;

vecWDVECIRCnew = sum([a1 a2 b1 b2 c2 c3].*matWCOEFF,2);
vecWDVECIRCnew(idx_flp,:) = vecWDVECIRCnew(idx_flp,:).*-1; 

end

