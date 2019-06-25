function [matVLST, matCENTER, matNEWWAKE, matCONTROL, matKINCON_P, matUINF, matVUINF, matPLEX, matDVECT, matROTANG] = fcnMOVEROTOR(matUINF, matVUINF, valRPM, valJ, valDIAM, valALPHA, valDELTIME, matVLST, matELST, vecTE, matDVE, matCENTER)

% Old trailing edge vertices
old_te = matVLST(matELST(vecTE,:),:);

%%
vecROTORRADPS = valRPM.*2.*pi./60;
vecROTORDEL = vecROTORRADPS.*valDELTIME;
dcmROTORSTEP = angle2dcm(vecROTORDEL,0,0,'ZXY');

translation = valJ.*(valRPM./60).*valDIAM.*fcnUINFWING(valALPHA, 0);

matVLST = matVLST*dcmROTORSTEP + translation.*valDELTIME;
matCENTER = matCENTER*dcmROTORSTEP + translation.*valDELTIME;
matUINF = cross(repmat([0,0,-vecROTORRADPS],length(matCENTER(:,1)),1),matCENTER) + translation;
matVUINF = cross(repmat([0,0,-vecROTORRADPS],length(matVLST(:,1)),1),matVLST) + translation;


%% Updating geometry
% matWCENTER
matCENTER = (matVLST(matDVE(:,1),:) + matVLST(matDVE(:,2),:) + matVLST(matDVE(:,3),:))./3;

% matWPLEX, matWDVECT, matWROTANG
P = permute(reshape(matVLST(matDVE(:,:)',:)', 3, 3, []), [2 1 3]);
DNORM = cross((matVLST(matDVE(:,3),:) - matVLST(matDVE(:,1),:)), (matVLST(matDVE(:,2),:) - matVLST(matDVE(:,1),:)), 2);
DNORM = DNORM./sqrt(sum(DNORM.^2, 2));
[matPLEX, matDVECT, matROTANG] = fcnTRITOLEX(P, DNORM, matCENTER);

matKINCON_P = matCENTER;
matCONTROL = matCENTER;

% New trailing edge vertices
new_te = matVLST(matELST(vecTE,:),:);

%% Output
matNEWWAKE(:,:,2) = [new_te(1:end/2,:); new_te((end/2)+1:end,:)];
matNEWWAKE(:,:,3) = [new_te((end/2)+1:end,:); old_te(1:end/2,:)];
matNEWWAKE(:,:,1) = [old_te(1:end/2,:); old_te((end/2)+1:end,:)];

end

