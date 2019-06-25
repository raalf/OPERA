function [matWELST, matWVLST, matWDVE, valWNELE, matWEIDX, matWPLEX, matWDVECT, matWCENTER, matWROTANG, matWVGRID] = ...
    fcnRELAX5(valDELTIME, valNELE, matCOEFF, matPLEX, valWNELE, matWCOEFF, matWDVE, matWVLST, matWPLEX, valWSIZE, ...
    matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST, matWEIDX, vecDVESYM, vecWDVESYM, vecWSYM, matWVGRID, vecWDVEFLIP)

% Not moving some vertices
oldest_edge = matWEIDX((1:valWSIZE) + valWSIZE,3); % Trailing edge of oldest wake row
dont_move = unique([matWELST(vecWLE,:); matWELST(vecWTE,:); matWELST(oldest_edge,:)]);
move = true(size(matWVLST,1),1);
move(dont_move) = false;

%% Getting velocities at wake vertices
tmp2 = zeros(size(matWVLST));
tmp2(move,:) = fcnSDVEVEL(matWVLST(move,:), valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM, [], 5e-3) + fcnSDVEVEL(matWVLST(move,:), valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, vecWDVESYM, [], 5e-3);
tmp2(matWELST(vecWSYM,:),2) = 0; % No y-component on symmetry line

%% Moving
% Moving vertices
matWVLST = matWVLST + tmp2.*valDELTIME;

%% Recreating wake point matrix, and regenerating wake HDVE parameters
% matWCENTER
matWCENTER = (matWVLST(matWDVE(:,1),:) + matWVLST(matWDVE(:,2),:) + matWVLST(matWDVE(:,3),:))./3;

% matWPLEX, matWDVECT, matWROTANG
P = permute(reshape(matWVLST(matWDVE(:,:)',:)', 3, 3, []), [2 1 3]);
DNORM = cross(matWVLST(matWDVE(:,2),:) - matWVLST(matWDVE(:,3),:), matWVLST(matWDVE(:,1),:) - matWVLST(matWDVE(:,3),:), 2);
DNORM = DNORM./sqrt(sum(DNORM.^2,2));
DNORM(vecWDVEFLIP,:) = DNORM(vecWDVEFLIP,:).*-1;
[matWPLEX, matWDVECT, matWROTANG] = fcnTRITOLEX(P, DNORM, matWCENTER);

end


