function [matVLST, matCENTER, matNEWWAKE, matUINF, matVUINF, matUINF_KK, matPLEX, matDVECT, matROTANG, matSPANDIR, matKINCON_P, matROTORHUB, matROTORAXIS] = ...
    fcnMOVESURFACES(valNELE, valUINF, valROTORS, valALPHA, valBETA, vecROTORRPM, valDELTIME, matROTORHUB, matVLST, matELST, vecTE, matDVE, ...
    vecDVEFLIP, matSPANDIR, matKINCON_P, vecDVEROTOR, vecKINCON_DVE, matROTORAXIS, valROTSTART, vecROTRATE, valTIMESTEP)

% Old trailing edge vertices
old_te = matVLST(matELST(vecTE,:),:);

%% Translational movement
translation = fcnUINFWING(valALPHA, valBETA).*valUINF.*valDELTIME;

matKINCON_P = matKINCON_P - translation;
matVLST = matVLST - translation;

%% Rotation of the rotors about the hub
for n = 1:valROTORS
    matROTORHUB(n,:) = matROTORHUB(n,:) - translation;

    idxDVEROTOR = vecDVEROTOR == n;
    
    tmpV = matVLST(matDVE(idxDVEROTOR,:),:);
    idx_k = ismember(vecKINCON_DVE,find(idxDVEROTOR));
    tmpK = matKINCON_P(idx_k,:);
    
    tmpV = tmpV - matROTORHUB(n,:);
    tmpK = tmpK - matROTORHUB(n,:);

    % transform rotor from hub plane to xy plane
    dcmXY2HUB = quat2dcm(fcnAXANG2QUAT(vrrotvec([0 0 1], matROTORAXIS(n,:))));
    tmpV = tmpV/dcmXY2HUB;
    tmpK = tmpK/dcmXY2HUB;
    tmpSD = matSPANDIR(idxDVEROTOR,:)/dcmXY2HUB;
    
    % Move the rotor in the hub plane
    tmpROTORRADPS = vecROTORRPM(n).*pi./30;
    vecROTORDEL = tmpROTORRADPS.*valDELTIME;
    dcmROTORSTEP = angle2dcm(vecROTORDEL,0,0,'ZXY');
    
    tmpV = tmpV*dcmROTORSTEP;
    tmpK = tmpK*dcmROTORSTEP;
    tmpSD = tmpSD*dcmROTORSTEP;
        
    if valTIMESTEP > valROTSTART    
        tmpDCM = angle2dcm(vecROTRATE(3)*0*valDELTIME, vecROTRATE(2)*valDELTIME, vecROTRATE(1)*valDELTIME);
        matROTORAXIS(n,:) = matROTORAXIS(n,:)*tmpDCM;
        dcmXY2HUB = quat2dcm(fcnAXANG2QUAT(vrrotvec([0 0 1], matROTORAXIS(n,:))));   
    end
    
    % Transform rotor back to hub
    tmpV = tmpV*dcmXY2HUB;
    tmpK = tmpK*dcmXY2HUB;
    tmpSD = tmpSD*dcmXY2HUB;
    
    tmpV = tmpV + matROTORHUB(n,:);
    tmpK = tmpK + matROTORHUB(n,:);
    
    matVLST(matDVE(idxDVEROTOR,:),:) = tmpV;
    matKINCON_P(idx_k,:) = tmpK;
    matSPANDIR(idxDVEROTOR,:) = tmpSD;    
end

%% Updating geometry
% matWCENTER
matCENTER = (matVLST(matDVE(:,1),:) + matVLST(matDVE(:,2),:) + matVLST(matDVE(:,3),:))./3;

% matWPLEX, matWDVECT, matWROTANG
P = permute(reshape(matVLST(matDVE(:,:)',:)', 3, 3, []), [2 1 3]);
DNORM = cross((matVLST(matDVE(:,3),:) - matVLST(matDVE(:,1),:)), (matVLST(matDVE(:,2),:) - matVLST(matDVE(:,1),:)), 2);
DNORM = DNORM./sqrt(sum(DNORM.^2, 2));
DNORM(vecDVEFLIP,:) = DNORM(vecDVEFLIP,:).*-1;
[matPLEX, matDVECT, matROTANG] = fcnTRITOLEX(P, DNORM, matCENTER);

% New trailing edge vertices
new_te = matVLST(matELST(vecTE,:),:);

%% Output
matNEWWAKE(:,:,2) = [new_te(1:end/2,:); new_te((end/2)+1:end,:)];
matNEWWAKE(:,:,3) = [new_te((end/2)+1:end,:); old_te(1:end/2,:)];
matNEWWAKE(:,:,1) = [old_te(1:end/2,:); old_te((end/2)+1:end,:)];

%% New speeds
[matUINF, matVUINF, matUINF_KK] = fcnUINF(valUINF, valALPHA, valBETA, valNELE, valROTORS, matVLST, matKINCON_P, ...
    vecROTORRPM, vecDVEROTOR, matCENTER, matDVE, matROTORHUB, matROTORAXIS, vecKINCON_DVE, valROTSTART, vecROTRATE, valTIMESTEP);

end

