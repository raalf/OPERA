function [matWCOEFF] = fcnDWAKE(type, valTIMESTEP, strATYPE, vecULS, valWNELE, vecWLE, vecWLEDVE, vecWTE, vecWTEDVE, matWEATT, matWELST, matWROTANG, matWCENTER, matWVLST, vecTE, vecTEDVE, matCOEFF, matCENTER, matROTANG, matWCOEFF, matWPLEX, vecWDVECIRC, vecWSYM, vecWSYMDVE, vecWDVESYM, matWVATT, matWDVE, vecWOTE, vecWOTEDVE)

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

vort_le = [fcnDVORT2(pts(:,:,1), vecWLEDVE, valWNELE, matWCENTER, matWROTANG, 'A');...
    fcnDVORT2(pts(:,:,2), vecWLEDVE, valWNELE, matWCENTER, matWROTANG, 'A');...
    fcnDVORT2(pts(:,:,3), vecWLEDVE, valWNELE, matWCENTER, matWROTANG, 'A')];
% vort_le = fcnDVORT2(pts(:,:,3), vecWLEDVE, valWNELE, matWCENTER, matWROTANG, 'A');
res_vort_le = zeros(size(vort_le,1),1);
  

%% Integrated circulation
circ_int = [];
res_circ_int = [];
vort_steady = [];
res_steady = [];
vort_te = [];
res_vort_te = [];

if strcmpi(type, 'UNSTEADY')   
    pts = [];
    pts(:,:,1) = matWVLST(matWELST(vecWOTE,1),:);
    pts(:,:,2) = matWVLST(matWELST(vecWOTE,2),:);
    pts(:,:,3) = (pts(:,:,1) + pts(:,:,2))./2;
    vort_te = [fcnDVORT2(pts(:,:,1), vecWOTEDVE, valWNELE, matWCENTER, matWROTANG, 'A');...
    fcnDVORT2(pts(:,:,2), vecWOTEDVE, valWNELE, matWCENTER, matWROTANG, 'A');...
    fcnDVORT2(pts(:,:,3), vecWOTEDVE, valWNELE, matWCENTER, matWROTANG, 'A')];
    % vort_te = fcnDVORT2(pts(:,:,3), vecWOTEDVE, valWNELE, matWCENTER, matWROTANG, 'A');
    res_vort_te = zeros(size(vort_te,1),1);
elseif strcmpi(type, 'STEADY')
    vort_steady = fcnDVORT1([1:valWNELE]', valWNELE, 'A');
    res_steady = zeros(size(vort_steady,1),1); 
end

%%
DW = [circ; vort; circ_le; vort_le; circ_int; vort_steady; vort_te];
RW = [res_circ; res_vort; res_circ_le; res_vort_le; res_circ_int; res_steady; res_vort_te];

matWCOEFF = DW\RW;
matWCOEFF = reshape(matWCOEFF,6,valWNELE,1)';

end

