function [matWELST, matWVLST, matWDVE, valWNELE, matWEIDX, matWPLEX, matWDVECT, matWVATT, matWCENTER, matWROTANG, matWAKEGEOM, matWVGRID] = ...
    fcnRELAX5(valTIMESTEP, valDELTIME, valNELE, matCOEFF, matPLEX, valWNELE, matWCOEFF, matWDVE, matWVLST, matWPLEX, valWSIZE, ...
    matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST, matWVATT, matWEIDX, vecDVESYM, vecWDVESYM, vecWSYM, matWVGRID, valPRESTEPS, vecWDVEFLIP)

% Not moving some vertices
oldest_edge = matWEIDX((1:valWSIZE) + valWSIZE,3); % Trailing edge of oldest wake row
dont_move = unique([matWELST(vecWLE,:); matWELST(vecWTE,:); matWELST(oldest_edge,:)]);
move = true(size(matWVLST,1),1);
move(dont_move) = false;

%% Getting velocities at wake vertices
tmp2 = zeros(size(matWVLST));
tmp2(move,:) = fcnSDVEVEL(matWVLST(move,:), valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM, false(valNELE,1)) + fcnSDVEVEL(matWVLST(move,:), valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, vecWDVESYM, false(valWNELE,1));
tmp2(matWELST(vecWSYM,:),2) = 0; % No y-component on symmetry line

%%
% forward = [matWVGRID(1,:); matWVGRID(1:end-1,:)];
% behind = [matWVGRID(2:end,:); matWVGRID(end,:)];
% tmp2 = (tmp2(forward,:) + 2.*tmp2 + tmp2(behind,:))./4;
% tmp2(dont_move,:) = tmp2(dont_move,:).*0;

%% Moving
% Moving vertices
matWVLST = matWVLST + tmp2.*valDELTIME;

%% Recreating wake point matrix, and regenerating wake HDVE parameters
matWAKEGEOM(:,:,1) = matWVLST(matWDVE(:,1,1),:);
matWAKEGEOM(:,:,2) = matWVLST(matWDVE(:,2,1),:);
matWAKEGEOM(:,:,3) = matWVLST(matWDVE(:,3,1),:);

% matWCENTER
matWCENTER = (matWVLST(matWDVE(:,1),:) + matWVLST(matWDVE(:,2),:) + matWVLST(matWDVE(:,3),:))./3;

% matWPLEX, matWDVECT, matWROTANG
P = permute(reshape(matWVLST(matWDVE(:,:)',:)', 3, 3, []), [2 1 3]);
DNORM = cross(matWVLST(matWDVE(:,2),:) - matWVLST(matWDVE(:,3),:), matWVLST(matWDVE(:,1),:) - matWVLST(matWDVE(:,3),:), 2);
DNORM = DNORM./sqrt(sum(DNORM.^2,2));
DNORM(vecWDVEFLIP,:) = DNORM(vecWDVEFLIP,:).*-1;
[matWPLEX, matWDVECT, matWROTANG] = fcnTRITOLEX(P, DNORM, matWCENTER, 'WAKE', []);

end


