function [matWELST, matWVLST, matWDVE, valWNELE, matWEIDX, matWPLEX, matWDVECT, matWVATT, matWCENTER, matWROTANG, matWAKEGEOM, matWVGRID] = ...
    fcnRELAX5(valTIMESTEP, valDELTIME, valNELE, matCOEFF, matPLEX, valWNELE, matWCOEFF, matWDVE, matWVLST, matWPLEX, valWSIZE, ...
    matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST, matWVATT, matWEIDX, vecDVESYM, vecWDVESYM, vecWSYM, matWVGRID, valPRESTEPS)

%% Getting velocities at wake vertices
tmp2 = fcnSDVEVEL(matWVLST, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM) + fcnSDVEVEL(matWVLST, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, vecWDVESYM);

%% Moving
% Not moving some vertices
tmp2(matWELST(vecWLE,:),:) = 0; % Trailing edge of wing
tmp2(matWELST(vecWTE,:),:) = 0; % Trailing edge of wing
tmp2(matWELST(vecWSYM,:),2) = 0; % No y-component on symmetry line
oldest_edge = matWEIDX((1:valWSIZE) + valWSIZE,3); % Trailing edge of oldest wake row
tmp2(matWELST(oldest_edge,:),:) = 0;

% Moving vertices
matWVLST = matWVLST + tmp2.*valDELTIME;

%% Recreating wake point matrix, and regenerating wake HDVE parameters
matWAKEGEOM(:,:,1) = matWVLST(matWDVE(:,1,1),:);
matWAKEGEOM(:,:,2) = matWVLST(matWDVE(:,2,1),:);
matWAKEGEOM(:,:,3) = matWVLST(matWDVE(:,3,1),:);

matWETA = nan(valWNELE, 1);
matWETA(reshape([1:valWSIZE.*2:valWNELE]' + [1:valWSIZE]-1, [], 1),1) = 1;
matWETA(reshape([1:valWSIZE.*2:valWNELE]' + [1:valWSIZE] + valWSIZE - 1, [], 1),1) = 2;

% matWCENTER
matWCENTER = (matWVLST(matWDVE(:,1),:) + matWVLST(matWDVE(:,2),:) + matWVLST(matWDVE(:,3),:))./3;

% matWPLEX, matWDVECT, matWROTANG
P = permute(reshape(matWVLST(matWDVE(:,:)',:)', 3, 3, []), [2 1 3]);
DNORM = cross((matWVLST(matWDVE(:,3),:) - matWVLST(matWDVE(:,1),:)), (matWVLST(matWDVE(:,2),:) - matWVLST(matWDVE(:,1),:)), 2);
DNORM = DNORM./sqrt(sum(DNORM.^2, 2));
[matWPLEX, matWDVECT, matWROTANG] = fcnTRITOLEX(P, DNORM, matWCENTER, 'WAKE', matWETA);



end


