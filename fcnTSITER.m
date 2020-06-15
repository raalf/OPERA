function [vecR, matCOEFF, vecVMU, vecEMU, matWCOEFF, vecWVMU, vecWEMU, vecTSITER, vecRSQUARED] = ...
    fcnTSITER(matPLEX, matCOEFF, valDLEN, valTIMESTEP, matUINF_KK, valWNELE, ...
    matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matKINCON_P, matKINCON_DVE, ...
    matDVECT, matD, vecR, valNELE, matVLST, matVATT, matCENTER, matROTANG, ...
    matDVE, matELST, matEIDX, strATYPE, vecWLE, matWVGRID, matWEGRID, matWE2GRID, ...
    vecWVMU, vecWEMU, matWELST, matWVLST, vecTEDVE, vecWOTE, matWDVE, matWEIDX, ...
    matEATT, boolKINCON, vecTE, vecVMU, vecEMU, valAZPERREV, valZTOL, vecRSQUARED, ...
    vecDVEFLIP, vecWDVEFLIP, vecHUB, vecWOFF, matWOFF, matWNORM)

count = 0;
delt = 1;
delt_old = 5;
sigma = 1; % relaxation

if any(~isnan(vecWOFF))
    for k = 1:length(vecWOFF)
        [tmpROTANG(:,:,k), tmpCENTER(:,:,k), tmpWROTANG(:,:,k), tmpWCENTER(:,:,k)] = ...
            fcnMIRROR(matVLST, matDVE, vecDVEFLIP, matWVLST, matWDVE, vecWDVEFLIP, ...
            vecHUB, vecWOFF(k), matWOFF(k,:), matWNORM(k,:));
    end
end

if valTIMESTEP > 1 && strcmpi(strATYPE{3},'UNSTEADY')
    idx = 1:(valWNELE - valWSIZE*2);
    w_ind_ow = (fcnSDVEVEL(matKINCON_P + matDVECT(matKINCON_DVE,:,3).*valZTOL, length(idx), matWCOEFF(idx,:), matWPLEX(:,:,idx), matWROTANG(idx,:), matWCENTER(idx,:), [], [], 0) + ...
        fcnSDVEVEL(matKINCON_P - matDVECT(matKINCON_DVE,:,3).*valZTOL, length(idx), matWCOEFF(idx,:), matWPLEX(:,:,idx), matWROTANG(idx,:), matWCENTER(idx,:), [], [], 0))./2;
    
    if any(~isnan(vecWOFF))
        w_ind_mirror = zeros(size(w_ind_ow));
        for k = 1:length(vecWOFF)
            w_ind_mirror = w_ind_mirror + ...
                fcnSDVEVEL(matKINCON_P, valNELE, matCOEFF, matPLEX, tmpROTANG(:,:,k), tmpCENTER(:,:,k), [], [], 0) + ...
                fcnSDVEVEL(matKINCON_P, length(idx), matWCOEFF(idx,:), matWPLEX(:,:,idx), tmpWROTANG(idx,:,k), tmpWCENTER(idx,:,k), [], [], 0);
        end
        w_ind_ow = w_ind_ow + w_ind_mirror;
    end
    
end

while ~isnan(delt) && delt >= 0.0001
    
    r2_1 = [vecVMU; vecEMU];
    
    if strcmpi(strATYPE{3},'UNSTEADY')
        idx = ((valWNELE - valWSIZE*2) + 1):valWNELE;
        w_ind_in = (fcnSDVEVEL(matKINCON_P + matDVECT(matKINCON_DVE,:,3).*valZTOL, length(idx), matWCOEFF(idx,:), matWPLEX(:,:,idx), matWROTANG(idx,:), matWCENTER(idx,:), [], [], 0) + ...
            fcnSDVEVEL(matKINCON_P - matDVECT(matKINCON_DVE,:,3).*valZTOL, length(idx), matWCOEFF(idx,:), matWPLEX(:,:,idx), matWROTANG(idx,:), matWCENTER(idx,:), [], [], 0))./2;
        
        if any(~isnan(vecWOFF))
            w_ind_mirror = zeros(size(w_ind_in));
            for k = 1:length(vecWOFF)
                w_ind_mirror = w_ind_mirror + fcnSDVEVEL(matKINCON_P, length(idx), matWCOEFF(idx,:), matWPLEX(:,:,idx), tmpWROTANG(idx,:,k), tmpWCENTER(idx,:,k), [], [], 0);
            end
            w_ind_in = w_ind_in + w_ind_mirror;
        end
        
        if valTIMESTEP > 1
            w_ind_in = w_ind_in + w_ind_ow;
        end
    else
        w_ind_in = [];
    end
    
    vecR = fcnRWING(strATYPE, valZTOL, valDLEN, valTIMESTEP, matUINF_KK, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, matKINCON_P, matKINCON_DVE, matDVECT, ...
        matVLST, matDVE, matCOEFF, matPLEX, vecDVEFLIP, matWVLST, matWDVE, vecWDVEFLIP, vecHUB, vecWOFF, matWOFF, matWNORM, w_ind_in);
    
    matCOEFF_new = fcnSOLVED(matD, vecR, valNELE);
    matCOEFF = sigma.*(matCOEFF_new - matCOEFF) + matCOEFF;
    
    [vecVMU, vecEMU] = fcnVEMU(matVLST, matVATT, matCENTER, matROTANG, matCOEFF, matELST, matEATT, vecTE);
    matCOEFF = fcnADJCOEFF(vecVMU, vecEMU, matVLST, matCENTER, matROTANG, matDVE, matCOEFF, matELST, matEIDX, valNELE);
    
    % Update wake coefficients
    [vecWVMU, vecWEMU] = fcnWAKEMU(strATYPE, vecWLE, matWVGRID, matWEGRID, matWE2GRID, vecWVMU, vecWEMU, matWELST, matWVLST, vecTEDVE, matCOEFF, matCENTER, matROTANG, vecWOTE, valAZPERREV, valTIMESTEP);
    matWCOEFF = fcnADJCOEFF(vecWVMU, vecWEMU, matWVLST, matWCENTER, matWROTANG, matWDVE, matWCOEFF, matWELST, matWEIDX, valWNELE);
    
    r2_2 = [vecVMU; vecEMU];
    delt = max(abs(r2_2 - r2_1));
        
    count = count + 1;
    disp(['       ',num2str(count),': ', num2str(delt)])
    
    if delt > delt_old
        break;
    end
    delt_old = delt;
end
vecTSITER(valTIMESTEP) = count;
[vecRSQUARED(valTIMESTEP,1), vecRSQUARED(valTIMESTEP,2)] = fcnRSQUARED(matD, vecR, matCOEFF, boolKINCON, matVLST, matVATT, matCENTER, matROTANG, matELST, matEATT);

end