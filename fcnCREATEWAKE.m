function [matWAKEGEOM, matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC,...
    matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWCOEFF, matWROTANG] = fcnCREATEWAKE(valTIMESTEP, matNEWWAKE, matWAKEGEOM, matCOEFF, valWSIZE, ...
    vecTE, vecTEDVE, matCENTER, matROTANG, matWCOEFF)

if valTIMESTEP <= 1
    matWAKEGEOM = matNEWWAKE;
else
    matWAKEGEOM = cat(1, matWAKEGEOM, matNEWWAKE);
end

[~, matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG] = fcnTRIANG(matWAKEGEOM);
vecWLEDVE = [(valWNELE - 2*valWSIZE + 1):(valWNELE - valWSIZE)]'; % Post trailing edge row of wake HDVEs
vecWLE = matWEIDX(vecWLEDVE,2);
vecWTEDVE = [(valWNELE - valWSIZE + 1):valWNELE]';
vecWTE = matWEIDX(vecWTEDVE,3);

if valTIMESTEP <= 1
    matWCOEFF = fcnDWAKENEW(valTIMESTEP, valWNELE, vecWLE, vecWLEDVE, vecWTE, vecWTEDVE, matWEATT, matWELST, matWROTANG, matWCENTER, matWVLST, vecTE, vecTEDVE, matCOEFF, matCENTER, matROTANG, []);
else
    tmp = fcnDWAKENEW(valTIMESTEP, valWNELE, vecWLE, vecWLEDVE, vecWTE, vecWTEDVE, matWEATT, matWELST, matWROTANG, matWCENTER, matWVLST, vecTE, vecTEDVE, matCOEFF, matCENTER, matROTANG, matWCOEFF);
    matWCOEFF = [matWCOEFF; tmp];
end




