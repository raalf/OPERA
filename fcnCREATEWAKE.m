function [matWAKEGEOM, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC,...
    matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWCOEFF, matWROTANG, vecWLE, vecWTE, ...
    vecWLEDVE, vecWTEDVE, vecWDVECIRC, vecWSYM, vecWSYMDVE, vecWDVESYM, matWVGRID, vecWOTE, ...
    vecWOTEDVE, matWEGRID, matWE2GRID, vecWVMU, vecWEMU] = fcnCREATEWAKE(valTIMESTEP, strATYPE, vecULS, matNEWWAKE, matWAKEGEOM, matCOEFF, valWSIZE, ...
    vecTE, vecTEDVE, matCENTER, matROTANG, matWCOEFF, matWPLEX, vecWDVECIRC, vecWSYMDVE, vecSYMDVE, vecWDVESYM, vecDVESYM, vecWSYM, matWVGRID, matWVLST, matWELST, matWEGRID, vecWVMU, vecWEMU, vecWDVEFLIP)

if valTIMESTEP <= 1
    matWAKEGEOM = matNEWWAKE;
else
    matWAKEGEOM = cat(1, matWAKEGEOM, matNEWWAKE);
end

valWNELE = valTIMESTEP.*valWSIZE.*2;

matWETA = nan(valWNELE, 1);
matWETA(reshape([1:valWSIZE.*2:valWNELE]' + [1:valWSIZE]-1, [], 1),1) = 1;
matWETA(reshape([1:valWSIZE.*2:valWNELE]' + [1:valWSIZE] + valWSIZE - 1, [], 1),1) = 2;

old_WVLST = matWVLST;
old_WELST = matWELST;
[~, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG] = fcnTRIANG(matWAKEGEOM, 'WAKE', matWETA, vecWDVEFLIP);

vecWLEDVE = [(valWNELE - 2*valWSIZE + 1):(valWNELE - valWSIZE)]'; % Post trailing edge row of wake HDVEs
vecWLE = matWEIDX(vecWLEDVE,2);
vecWTEDVE = [(valWNELE - valWSIZE + 1):valWNELE]';
vecWTE = matWEIDX(vecWTEDVE,3);

vecWOTEDVE = [(1:valWSIZE) + valWSIZE]';
vecWOTE = matWEIDX(vecWOTEDVE, 3);

if any(vecDVESYM)
    vecWDVESYM = [vecWDVESYM; repmat(vecDVESYM(vecTEDVE), 2, 1)];
end
idx = ismember(vecTEDVE, vecSYMDVE);
vecWSYMDVE = [vecWSYMDVE; vecWLEDVE(idx)];
vecWSYM = [vecWSYM; matWEIDX(vecWSYMDVE, 1)];

%% Wake circulation grid points
% fcnPLOTWAKE(1, gcf, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, valWSIZE, 0);
% view([90 90])
% hFig1 = fcnPLOTWAKE(1, gcf, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, valWSIZE, 0, matWVGRID);

if valTIMESTEP == 1
    matWVGRID = [unique(matWELST(vecWLE,:)', 'stable')'; unique(matWELST(vecWTE,:)', 'stable')'];
    matWEGRID = [vecWLE matWEIDX(vecWLEDVE,3) vecWTE]';
    vecWVMU = zeros(size(matWVLST,1),1);
    vecWEMU = zeros(size(matWELST,1),1);
else
    % find matWVGRID vertices from old_WVLST in new matWVLST
    [~,vmap] = ismember(old_WVLST, matWVLST, 'rows');
    [~,emap] = ismember(sort(vmap(old_WELST),2), sort(matWELST,2), 'rows');
    matWVGRID = [unique(matWELST(vecWLE,:)', 'stable')'; vmap(matWVGRID)];
    matWEGRID = [vecWLE'; matWEIDX(vecWLEDVE,3)'; emap(matWEGRID)];
    
    old_WVMU = vecWVMU;
    old_WEMU = vecWEMU;
    vecWVMU = zeros(size(matWVLST,1),1);
    vecWEMU = zeros(size(matWELST,1),1);
    
    if strcmpi(strATYPE{3},'UNSTEADY')
        vecWVMU(vmap) = old_WVMU;
        vecWEMU(emap) = old_WEMU;
    end
end

tmp = permute(cat(3,matWVGRID(1:end-1,:), matWVGRID(2:end,:)), [1 3 2]);
tmp = reshape(permute(tmp, [2 1 3]), size(tmp, 2), [])';
[~,idx] = ismember(sort(tmp,2), sort(matWELST,2),'rows');
matWE2GRID = reshape(idx, valTIMESTEP, size(matWVGRID,2));

%% Coefficients
[vecWVMU, vecWEMU] = fcnWAKEMU(strATYPE, vecWLE, matWVGRID, matWEGRID, matWE2GRID, vecWVMU, vecWEMU, matWELST, matWVLST, vecTEDVE, matCOEFF, matCENTER, matROTANG);
matWCOEFF = fcnADJCOEFF(vecWVMU, vecWEMU, matWVLST, matWCENTER, matWROTANG, matWDVE, matWCOEFF, matWELST, matWEIDX, valWNELE);


end




