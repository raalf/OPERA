function [matVLST, matCENTER, matNEWWAKE, matCONTROL, matKINCON_P, matUINF, matVUINF, matPLEX, matDVECT, matROTANG, vecWDVEFLIP] = fcnMOVEROTOR(matUINF, matVUINF, valRPM, valJ, valDIAM, valALPHA, valDELTIME, matVLST, matELST, vecTE, matDVE, vecWDVEFLIP)

% Old trailing edge vertices
old_te = matVLST(matELST(vecTE,:),:);

matVLST = matVLST + matVUINF.*valDELTIME;

% matWCENTER
matCENTER = (matVLST(matDVE(:,1),:) + matVLST(matDVE(:,2),:) + matVLST(matDVE(:,3),:))./3;

% matWPLEX, matWDVECT, matWROTANG
P = permute(reshape(matVLST(matDVE(:,:)',:)', 3, 3, []), [2 1 3]);
DNORM = cross((matVLST(matDVE(:,3),:) - matVLST(matDVE(:,1),:)), (matVLST(matDVE(:,2),:) - matVLST(matDVE(:,1),:)), 2);
DNORM = DNORM./sqrt(sum(DNORM.^2, 2));
[matPLEX, matDVECT, matROTANG] = fcnTRITOLEX(P, DNORM, matCENTER, 'SURFACE');

matKINCON_P = matCENTER;
matCONTROL = matCENTER;

% Update matVUINF, matUINF
rps = valRPM./60;
rev = 2*pi*rps*valDELTIME;

dist = 2.*pi.*(sqrt(sum(matCENTER(:,1:2).^2,2))).*(rev/(2*pi));
dir = cross(repmat([0 0 1], size(matCENTER,1), 1), matCENTER, 2);
dir = dir./sqrt(sum(dir.^2, 2));
matUINF = dir.*(dist./valDELTIME) + valJ.*rps.*valDIAM.*fcnUINFWING(valALPHA, 0);

dist = 2.*pi.*(sqrt(sum(matVLST(:,1:2).^2,2))).*(rev/(2*pi));
dir = cross(repmat([0 0 1], size(matVLST,1), 1), matVLST, 2);
dir = dir./sqrt(sum(dir.^2, 2));
matVUINF = dir.*(dist./valDELTIME) + valJ.*rps.*valDIAM.*fcnUINFWING(valALPHA, 0);

% New trailing edge vertices
new_te = matVLST(matELST(vecTE,:),:);

matNEWWAKE(:,:,2) = [new_te(1:end/2,:); new_te((end/2)+1:end,:)];
matNEWWAKE(:,:,3) = [new_te((end/2)+1:end,:); old_te(1:end/2,:)];
matNEWWAKE(:,:,1) = [old_te(1:end/2,:); old_te((end/2)+1:end,:)];

end

