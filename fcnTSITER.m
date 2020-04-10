function [vecR, matCOEFF, vecVMU, vecEMU, matWCOEFF, vecWVMU, vecWEMU, vecTSITER, vecRSQUARED] = ...
    fcnTSITER(matPLEX, matCOEFF, matDVEGRID, valDLEN, valTIMESTEP, matUINF_KK, valWNELE, ...
    matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matKINCON_P, matKINCON_DVE, ...
    matDVECT, vecWDVESYM, matD, vecR, valNELE, matVLST, matVATT, matCENTER, matROTANG, ...
    matDVE, matELST, matEIDX, strATYPE, vecWLE, matWVGRID, matWEGRID, matWE2GRID, ...
    vecWVMU, vecWEMU, matWELST, matWVLST, vecTEDVE, vecWOTE, matWDVE, matWEIDX, ...
    matEATT, boolKINCON, vecTE)


count = 0;
delt = 1;
sigma = 0.5; % relaxation
while ~isnan(delt) && delt >= 0.001
    
    [~, r2_1] = fcnRSQUARED(matD, vecR, matCOEFF, boolKINCON, matVLST, matVATT, matCENTER, matROTANG, matELST, matEATT);
        
    [vecR, wind] = fcnRWING(valDLEN, valTIMESTEP, matUINF_KK, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matKINCON_P, matKINCON_DVE, matDVECT, vecWDVESYM);
    matCOEFF_new = fcnSOLVED(matD, vecR, valNELE);
    matCOEFF = sigma.*(matCOEFF_new - matCOEFF) + matCOEFF;
    
    [vecVMU, vecEMU] = fcnVEMU(matVLST, matVATT, matCENTER, matROTANG, matCOEFF, matELST, matEATT, vecTE);
    matCOEFF = fcnADJCOEFF(vecVMU, vecEMU, matVLST, matCENTER, matROTANG, matDVE, matCOEFF, matELST, matEIDX, valNELE);
    
    % Update wake coefficients
    [vecWVMU, vecWEMU] = fcnWAKEMU(strATYPE, vecWLE, matWVGRID, matWEGRID, matWE2GRID, vecWVMU, vecWEMU, matWELST, matWVLST, vecTEDVE, matCOEFF, matCENTER, matROTANG, vecWOTE);
    matWCOEFF = fcnADJCOEFF(vecWVMU, vecWEMU, matWVLST, matWCENTER, matWROTANG, matWDVE, matWCOEFF, matWELST, matWEIDX, valWNELE);
        
    [~, r2_2] = fcnRSQUARED(matD, vecR, matCOEFF, boolKINCON, matVLST, matVATT, matCENTER, matROTANG, matELST, matEATT);
    delt = abs((r2_2 - r2_1)./r2_1);
    
    count = count + 1;
    
    if count > 50
        disp(['          Diverged. Error: ', num2str(delt), ' Count = ', num2str(count)]);
        break;
    end
    disp(['       ',num2str(count),': ', num2str(delt)])
end
vecTSITER(valTIMESTEP) = count;
[vecRSQUARED(valTIMESTEP,1), vecRSQUARED(valTIMESTEP,2)] = fcnRSQUARED(matD, vecR, matCOEFF, boolKINCON, matVLST, matVATT, matCENTER, matROTANG, matELST, matEATT);

end