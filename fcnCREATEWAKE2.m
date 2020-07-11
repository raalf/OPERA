function [matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWPLEX, ...
    matWDVECT, matWCENTER, matWCOEFF, matWROTANG, vecWLE, vecWTE, vecWLEDVE, vecWTEDVE, ...
    matWVGRID, vecWOTE, vecWOTEDVE, matWEGRID, matWE2GRID, vecWVMU, vecWEMU, vecWDVEFLIP, ...
    vecWDVESURFACE, matWDVEGRID] = fcnCREATEWAKE2(valTIMESTEP, flgSTEADY, matNEWWAKE, matCOEFF, valWSIZE, ...
    vecTEDVE, matCENTER, matROTANG, matWCOEFF, matWPLEX, matWVGRID, matWVLST, matWELST, matWEGRID, ...
    vecWVMU, vecWEMU, vecWDVEFLIP, matWCENTER, matWROTANG, matWDVECT, ...
    matWDVE, matWEIDX, vecWLEDVE, matWEATT, vecDVESURFACE, vecWDVESURFACE, matWDVEGRID)

valWNELE = valTIMESTEP*valWSIZE*2;
WDVEFLIP = [false(size(matNEWWAKE,1)./2, 1); true(size(matNEWWAKE,1)./2, 1)];
[~, WELST, WVLST, WDVE, ~, WEATT, WEIDX, WPLEX, WDVECT, ~, ~, WCENTER, WROTANG] = fcnTRIANG(matNEWWAKE, WDVEFLIP);

num_dves = size(matWDVE,1);
if valTIMESTEP > 1
    [idx1,b1] = ismember(WVLST, matWVLST, 'rows');
    
    % WDVE, WELST
    idx_WDVE = false(size(WDVE));
    idx_WELST = false(size(WELST));
    new = find(idx1);
    old = b1(idx1);
    for j = 1:length(new)
        tmp_DVE = WDVE == new(j) & ~idx_WDVE;
        WDVE(tmp_DVE) = old(j);
        idx_WDVE = idx_WDVE | tmp_DVE;
        
        tmp_WELST = WELST == new(j) & ~idx_WELST;
        WELST(tmp_WELST) = old(j);
        idx_WELST = idx_WELST | tmp_WELST;
    end
    
    old = unique(WDVE(~idx_WDVE));
    new_vt = [(size(matWVLST,1)+1):size(matWVLST,1) + length(new)];
    for j = 1:length(old)
        tmp_DVE = WDVE == old(j) & ~idx_WDVE;
        WDVE(tmp_DVE) = new_vt(j);
        
        tmp_WELST = WELST == old(j) & ~idx_WELST;
        WELST(tmp_WELST) = new_vt(j);
    end
    
    matWDVE = cat(1, matWDVE, WDVE);
    matWVLST = cat(1, matWVLST, WVLST(~idx1, :));
    
    % WEIDX
    [idx2,b2] = ismember(WELST, matWELST, 'rows');
    idx_WEIDX = false(size(WEIDX));
    new = find(idx2);
    old = b2(idx2);
    for j = 1:length(new)
        tmp_WEIDX = WEIDX == new(j);
        WEIDX(tmp_WEIDX) = old(j);
        idx_WEIDX = idx_WEIDX | tmp_WEIDX;
    end
    
    old = unique(WEIDX(~idx_WEIDX));
    new = [(size(matWELST,1)+1):size(matWELST,1) + length(old)]';
    for j = 1:length(old)
        tmp_WEIDX = WEIDX == old(j) & ~idx_WEIDX;
        WEIDX(tmp_WEIDX) = new(j);
    end
    
    matWEIDX = cat(1, matWEIDX, WEIDX);
    matWELST = cat(1, matWELST, WELST(~idx2, :));
else
    matWDVE = WDVE;
    matWVLST = WVLST;
    matWELST = WELST;
    matWEIDX = WEIDX;
    matWEATT = WEATT;
end

matWDVEGRID = [[1:valWSIZE, (valWSIZE+1):valWSIZE*2] + valWSIZE.*2.*(valTIMESTEP-1); matWDVEGRID];
matWCENTER = cat(1, matWCENTER, WCENTER);
matWROTANG = cat(1, matWROTANG, WROTANG);
matWDVECT = cat(1, matWDVECT, WDVECT);
matWPLEX = cat(3, matWPLEX, WPLEX);
vecWDVEFLIP = cat(1, vecWDVEFLIP, WDVEFLIP);
vecWDVESURFACE = [vecWDVESURFACE; repmat(vecDVESURFACE(vecTEDVE), 2, 1)];

% WEATT needs special care when valTIMESTEP > 1

if valTIMESTEP > 1
    % remove duplicate edges (TE of wake row)
    WEATT(idx2,:) = [];
    WEATT(WEATT > 0) = WEATT(WEATT > 0) + num_dves;
    % rename new vertices to old
    % concatenate?
    matWEATT = cat(1, matWEATT, WEATT);
end

% Wake leading and trailing edge DVE identifications
if valTIMESTEP > 1
    old_le_dve = vecWLEDVE;
end
vecWLEDVE = [(valWNELE - 2*valWSIZE + 1):(valWNELE - valWSIZE)]'; % Post trailing edge row of wake HDVEs
vecWTEDVE = [(valWNELE - valWSIZE + 1):valWNELE]';

if valTIMESTEP > 1
    matWEATT(matWEIDX(old_le_dve,2),:) = [old_le_dve vecWLEDVE];
end

% Wake leading and trailing edge DVE identifications
vecWLE = matWEIDX(vecWLEDVE,2);
vecWTE = matWEIDX(vecWTEDVE,3);
vecWOTEDVE = [(1:valWSIZE) + valWSIZE]';
vecWOTE = matWEIDX(vecWOTEDVE, 3);

%% Wake circulation grid points
WVMU = zeros(size(matWVLST,1),1);
WEMU = zeros(size(matWELST,1),1);

if valTIMESTEP == 1
    matWVGRID = [unique(reshape(matWDVE(vecWLEDVE,2:3)',1,[],1)','stable')'; unique(reshape(matWDVE(vecWTEDVE,3:-2:1)',1,[],1)','stable')'];
    matWEGRID = [vecWLE matWEIDX(vecWLEDVE,3) vecWTE]';
    vecWVMU = WVMU;
    vecWEMU = WEMU;
else
    matWVGRID = [unique(reshape(matWDVE(vecWLEDVE,2:3)',1,[],1)','stable')'; matWVGRID];
    matWEGRID = [[vecWLE matWEIDX(vecWLEDVE,3)]'; matWEGRID];
    
    old_WVMU = vecWVMU;
    old_WEMU = vecWEMU;
    vecWVMU = zeros(size(matWVLST,1),1);
    vecWEMU = zeros(size(matWELST,1),1);
    
    if flgSTEADY == false
        vecWVMU(1:length(old_WVMU)) = old_WVMU;
        vecWEMU(1:length(old_WEMU)) = old_WEMU;
    end
end

tmp = permute(cat(3, matWVGRID(1:end-1,:), matWVGRID(2:end,:)), [1 3 2]);
tmp = reshape(permute(tmp, [2 1 3]), size(tmp, 2), [])';
[~,idx] = ismember(sort(tmp,2), sort(matWELST,2),'rows');
matWE2GRID = reshape(idx, valTIMESTEP, size(matWVGRID,2));

%% Coefficients
[vecWVMU, vecWEMU] = fcnWAKEMU(flgSTEADY, vecWLE, matWVGRID, matWEGRID, matWE2GRID, vecWVMU, vecWEMU, matWELST, matWVLST, vecTEDVE, matCOEFF, matCENTER, matROTANG, vecWOTE);
matWCOEFF = fcnADJCOEFF(vecWVMU, vecWEMU, matWVLST, matWCENTER, matWROTANG, matWDVE, matWCOEFF, matWELST, matWEIDX, valWNELE);


end




