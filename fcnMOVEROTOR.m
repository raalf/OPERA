function [matVLST, matCENTER, matNEWWAKE, matCONTROL, matKINCON_P, matUINF, matVUINF, matPLEX, matDVECT, matROTANG] = fcnMOVEROTOR(matUINF, matVUINF, valRPM, valJ, valDIAM, valALPHA, valDELTIME, matVLST, matELST, vecTE, matDVE)

% 3D matrix of TE points, 3x?x2, rows are (x,y,z), columns are points, and depth is first and second TE point
te_flip(:,:,1) = matVLST(matELST(vecTE,1),:)';
te_flip(:,:,2) = matVLST(matELST(vecTE,2),:)';

idx_flip = [sqrt(sum(te_flip(1:2,:,1).^2,1)) > sqrt(sum(te_flip(1:2,:,2).^2,1))]'; % Finding out which trailing edges go from out to in, so we know what to flip
idx_flip2 = logical([zeros(length(idx_flip),1); idx_flip]); % idx_flip1 is the first point, idx_flip2 is the second point

% Old trailing edge vertices
old_te = matVLST(matELST(vecTE,:),:);

% Flipping the necessary vertices
temp = old_te(idx_flip,:);
old_te(idx_flip,:) = old_te(idx_flip2,:);
old_te(idx_flip2,:) = temp;

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

% Flipping the necessary vertices
temp = new_te(idx_flip,:);
new_te(idx_flip,:) = new_te(idx_flip2,:);
new_te(idx_flip2,:) = temp;

% Now that everything is flipped, the wake SHOULD generate with eta in the streamwise direction and
% the normals pointed upwards. If it doesn't, then this needs to be rewritten to ensure that is the case

% These vertices will be used to calculate the wake HDVE geometry
matNEWWAKE(:,:,2) = [new_te(1:end/2,:); new_te((end/2)+1:end,:)];
matNEWWAKE(:,:,3) = [new_te((end/2)+1:end,:); old_te((end/2)+1:end,:)];
matNEWWAKE(:,:,1) = [old_te(1:end/2,:); old_te(1:end/2,:)];

end

