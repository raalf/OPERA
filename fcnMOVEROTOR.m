function [matVLST, matCENTER, matNEWWAKE, matUINF, matVUINF, matPLEX, matDVECT, matROTANG, vecHUB, matSPANDIR] = fcnMOVEROTOR(atype, valRPM, valJ, valDIAM, valALPHA, valDELTIME, matVLST, matELST, vecTE, matDVE, matCENTER, vecDVEFLIP, vecHUB, matSPANDIR)

% Old trailing edge vertices
old_te = matVLST(matELST(vecTE,:),:);

%%
vecROTORRADPS = valRPM.*pi./30;
vecROTORDEL = vecROTORRADPS.*valDELTIME;
dcmROTORSTEP = angle2dcm(vecROTORDEL,0,0,'ZXY');

if strcmpi(atype, 'ROTOR')
    translation = valJ.*(valRPM.*(pi/30)).*(valDIAM/2).*fcnUINFWING(valALPHA, 0);
elseif strcmpi(atype, 'PROPELLER')
    translation = valJ.*(valRPM/60).*(valDIAM).*fcnUINFWING(valALPHA, 0);
end
    
tmpVLST = matVLST - vecHUB;
tmpVLST = tmpVLST*dcmROTORSTEP;

matSPANDIR = matSPANDIR*dcmROTORSTEP;

matUINF = cross(repmat([0,0,-vecROTORRADPS],length(matCENTER(:,1)),1),matCENTER - vecHUB) - translation;
matVUINF = cross(repmat([0,0,-vecROTORRADPS],length(matVLST(:,1)),1),matVLST - vecHUB) - translation;

vecHUB = vecHUB + translation.*valDELTIME;
matVLST = tmpVLST + vecHUB;

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

end

