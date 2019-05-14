function [matWAKEGEOM, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC,...
    matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWCOEFF, matWROTANG, vecWLE, vecWTE, vecWLEDVE, vecWTEDVE, vecWDVECIRC, vecWSYM, vecWSYMDVE, vecWDVESYM, matWVGRID] = fcnCREATEWAKE(valTIMESTEP, strATYPE, vecULS, matNEWWAKE, matWAKEGEOM, matCOEFF, valWSIZE, ...
    vecTE, vecTEDVE, matCENTER, matROTANG, matWCOEFF, matWPLEX, vecWDVECIRC, vecWSYMDVE, vecSYMDVE, vecWDVESYM, vecDVESYM, vecWSYM)

if valTIMESTEP <= 1
    matWAKEGEOM = matNEWWAKE;
else
    matWAKEGEOM = cat(1, matWAKEGEOM, matNEWWAKE);
end

valWNELE = valTIMESTEP.*valWSIZE.*2;

matWETA = nan(valWNELE, 1);
matWETA(reshape([1:valWSIZE.*2:valWNELE]' + [1:valWSIZE]-1, [], 1),1) = 1;
matWETA(reshape([1:valWSIZE.*2:valWNELE]' + [1:valWSIZE] + valWSIZE - 1, [], 1),1) = 2;

[~, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG] = fcnTRIANG(matWAKEGEOM, 'WAKE', matWETA);
vecWLEDVE = [(valWNELE - 2*valWSIZE + 1):(valWNELE - valWSIZE)]'; % Post trailing edge row of wake HDVEs
vecWLE = matWEIDX(vecWLEDVE,2);
vecWTEDVE = [(valWNELE - valWSIZE + 1):valWNELE]';
vecWTE = matWEIDX(vecWTEDVE,3);

if any(vecDVESYM)
    vecWDVESYM = [vecWDVESYM; repmat(vecDVESYM(vecTEDVE), 2, 1)];
end
idx = ismember(vecTEDVE, vecSYMDVE);
vecWSYMDVE = [vecWSYMDVE; vecWLEDVE(idx)];
vecWSYM = [vecWSYM; matWEIDX(vecWSYMDVE, 1)];

if valTIMESTEP <= 1
    [matWCOEFF, vecWDVECIRC] = fcnDWAKENEW(valTIMESTEP, strATYPE, vecULS, valWNELE, vecWLE, vecWLEDVE, vecWTE, vecWTEDVE, matWEATT, matWELST, matWROTANG, matWCENTER, matWVLST, vecTE, vecTEDVE, matCOEFF, matCENTER, matROTANG, [], matWPLEX, vecWSYM, vecWSYMDVE, vecWDVESYM, matWVATT, matWDVE);
%     matWCOEFF = fcnADJCOEFF(matWVLST, matWVATT, matWCENTER, matWROTANG, matWDVE, matWCOEFF, matWELST, matWEATT, matWEIDX, [], valWNELE, []);
else
    [tmp1, tmp2] = fcnDWAKENEW(valTIMESTEP, strATYPE, vecULS, valWNELE, vecWLE, vecWLEDVE, vecWTE, vecWTEDVE, matWEATT, matWELST, matWROTANG, matWCENTER, matWVLST, vecTE, vecTEDVE, matCOEFF, matCENTER, matROTANG, matWCOEFF, matWPLEX, vecWSYM, vecWSYMDVE, vecWDVESYM, matWVATT, matWDVE);
    matWCOEFF = [matWCOEFF; tmp1];
    vecWDVECIRC = [vecWDVECIRC; tmp2];
%     matWCOEFF = fcnADJCOEFF(matWVLST, matWVATT, matWCENTER, matWROTANG, matWDVE, matWCOEFF, matWELST, matWEATT, matWEIDX, [], valWNELE, []);
end

dvegrid = flipud(repmat([1:valWSIZE, valWSIZE*2], valTIMESTEP, 1) + [0:(valWSIZE*2):(valWSIZE*valTIMESTEP*2 - 1)]');
matWVGRID = [reshape(matWDVE(dvegrid,2), size(dvegrid)); reshape(matWDVE(dvegrid(end,1:end-1),1), size(dvegrid(end,1:end-1))) matWDVE(dvegrid(end,end), 3)];


end




