function [matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC,...
    matWPLEX, matWDVECT, matWCENTER, matWCOEFF, matWROTANG, vecWLE, vecWTE, ...
    vecWLEDVE, vecWTEDVE, vecWDVECIRC, vecWSYM, vecWSYMDVE, vecWDVESYM, matWVGRID, vecWOTE, ...
    vecWOTEDVE, matWEGRID, matWE2GRID, vecWVMU, vecWEMU, vecWDVEFLIP] = fcnCREATEWAKE2(valTIMESTEP, strATYPE, matNEWWAKE, matCOEFF, valWSIZE, ...
    vecTEDVE, matCENTER, matROTANG, matWCOEFF, matWPLEX, vecWDVECIRC, vecWSYMDVE, vecSYMDVE, vecWDVESYM, vecDVESYM, ...
    vecWSYM, matWVGRID, matWVLST, matWELST, matWEGRID, vecWVMU, vecWEMU, vecWDVEFLIP, matWCENTER, matWROTANG, matWDVECT, ...
    matWELOC, matWDVE, matWEIDX, vecWLEDVE, vecWTEDVE, matWEATT)

WDVEFLIP = [false(size(matNEWWAKE,1)./2, 1); true(size(matNEWWAKE,1)./2, 1)];

valWNELE = valTIMESTEP*valWSIZE*2;
[~, WELST, WVLST, WDVE, ~, WEATT, WEIDX, WELOC, WPLEX, WDVECT, ~, ~, WCENTER, WROTANG] = fcnTRIANG(matNEWWAKE, WDVEFLIP);

% Requiring no modification, just concatination
if valTIMESTEP > 1
    [idx1,~] = ismember(WVLST, matWVLST, 'rows');
    newWVLST = [WVLST(~idx1,:); WVLST(idx1,:)];
    [~,vmap] = ismember(WVLST, newWVLST, 'rows');
    
    WDVE = vmap(WDVE);
    WELST = vmap(WELST);

    [idx1,b] = ismember(matWVLST, newWVLST, 'rows');
    idx3 = find(~idx1);

    
    [idx2,b] = ismember(WDVE, idx3);
    WDVE = WDVE + size(matWVLST,1)/2;
    WDVE(b > 0) = b(b > 0);
    
    [~,b] = ismember(WELST, idx3);
    WELST(all(b,2),:) = [];
    b(all(b,2),:) = [];
    WELST = WELST + size(matWVLST,1);
    WELST(b > 0) = b(b > 0);
    
    [~,b] = ismember(WEIDX, idx3);
    WEIDX = WEIDX + size(matWVLST,1);
    WEIDX(b > 0) = b(b > 0);
    
    matWVLST = cat(1, matWVLST, WVLST(~idx1, :));
else
    matWVLST = WVLST;
end

hold on
fcnPLOTWAKE(1, gcf, WDVE, valWSIZE*2, matWVLST, WELST, WDVECT, WCENTER, valWSIZE, 0);
view([90 90])
hold off

matWCENTER = cat(1, matWCENTER, WCENTER);
matWROTANG = cat(1, matWROTANG, WROTANG);
matWDVECT = cat(1, matWDVECT, WDVECT);
matWPLEX = cat(3, matWPLEX, WPLEX);
matWELOC = cat(1, matWELOC, WELOC);
vecWDVEFLIP = cat(1, vecWDVEFLIP, WDVEFLIP);
matWDVE = cat(1, matWDVE, WDVE);
matWELST = cat(1, matWELST, WELST);
matWEIDX = cat(1, matWEIDX, WEIDX);

% WEATT needs special care when valTIMESTEP > 1
matWEATT = cat(1, matWEATT, WEATT);

% Wake leading and trailing edge DVE identifications
if valTIMESTEP > 1
    old_le_dve = vecWLEDVE;
end
vecWLEDVE = [(valWNELE - 2*valWSIZE + 1):(valWNELE - valWSIZE)]'; % Post trailing edge row of wake HDVEs
vecWTEDVE = [(valWNELE - valWSIZE + 1):valWNELE]';

if valTIMESTEP > 1
    matWEATT(matWEIDX(vecWLEDVE,3),:) = sort(matWEATT(matWEIDX(vecWLEDVE,3),:), 2);
    matWEATT(matWEIDX(vecWLEDVE,1),1) = old_le_dve;
end

% Wake leading and trailing edge DVE identifications
vecWLE = matWEIDX(vecWLEDVE,2);
vecWTE = matWEIDX(vecWTEDVE,3);
vecWOTEDVE = [(1:valWSIZE) + valWSIZE]';
vecWOTE = matWEIDX(vecWOTEDVE, 3);

% Symmetry
if any(vecDVESYM)
    vecWDVESYM = [vecWDVESYM; repmat(vecDVESYM(vecTEDVE), 2, 1)];
end
idx = ismember(vecTEDVE, vecSYMDVE);
vecWSYMDVE = [vecWSYMDVE; vecWLEDVE(idx)];
vecWSYM = [vecWSYM; matWEIDX(vecWSYMDVE, 1)];

%% Wake circulation grid points
% fcnPLOTWAKE(1, gcf, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, valWSIZE, 0);
% view([90 90])

WVGRID = [unique(matWELST(vecWLE,:)', 'stable')'; unique(matWELST(vecWTE,:)', 'stable')'];
WEGRID = [vecWLE matWEIDX(vecWLEDVE,3) vecWTE]';
WVMU = zeros(size(matWVLST,1),1);
WEMU = zeros(size(matWELST,1),1);

if valTIMESTEP == 1
    matWVGRID = WVGRID;
    matWEGRID = WEGRID;
    vecWVMU = WVMU;
    vecWEMU = WEMU;
else
    matWVGRID = cat(1, matWVGRID, WVGRID);
    matWEGRID = cat(1, matWEGRID, WEGRID);
    
    old_WVMU = vecWVMU;
    old_WEMU = vecWEMU;
    vecWVMU = zeros(size(matWVLST,1),1);
    vecWEMU = zeros(size(matWELST,1),1);
    
    if strcmpi(strATYPE{3},'UNSTEADY')
        vecWVMU(vmap) = old_WVMU;
        vecWEMU(emap) = old_WEMU;
    end
end

tmp = permute(cat(3, matWVGRID(1:end-1,:), matWVGRID(2:end,:)), [1 3 2]);
tmp = reshape(permute(tmp, [2 1 3]), size(tmp, 2), [])';
[~,idx] = ismember(sort(tmp,2), sort(matWELST,2),'rows');
matWE2GRID = reshape(idx, valTIMESTEP, size(matWVGRID,2));

%% Coefficients
[vecWVMU, vecWEMU] = fcnWAKEMU(strATYPE, vecWLE, matWVGRID, matWEGRID, matWE2GRID, vecWVMU, vecWEMU, matWELST, matWVLST, vecTEDVE, matCOEFF, matCENTER, matROTANG);
matWCOEFF = fcnADJCOEFF(vecWVMU, vecWEMU, matWVLST, matWCENTER, matWROTANG, matWDVE, matWCOEFF, matWELST, matWEIDX, valWNELE);


end




