function [matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG, matWAKEGEOM] = ...
                fcnRELAX(valTIMESTEP, matUINF, valDELTIME, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, ...
                matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST, matWVATT, matWEIDX, vecDVESYM, vecWDVESYM, vecWSYM)

% Finding induced velocities at all wake vertices
q_ind = fcnSDVEVEL(matWCENTER, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER) + fcnSDVEVEL(matWCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER);

%% Moving wake vertices
% matWCENTER = matWCENTER + q_ind.*valDELTIME;
for i = 1:size(matWVATT,1)
   tmp2(i,:) = mean(q_ind(~isnan(matWVATT(i,:)),:),1);
end

tmp2(matWELST(vecWLE,:),:) = tmp2(matWELST(vecWLE,:),:).*0; % Trailing edge of wing
tmp2(matWELST(vecWTE,:),:) = tmp2(matWELST(vecWTE,:),:).*0; % Trailing edge of LE wake element row
oldest_edge = matWEIDX((1:valWSIZE) + valWSIZE,3); % Trailing edge of oldest wake row
tmp2(matWELST(oldest_edge,:),:) = tmp2(matWELST(oldest_edge,:),:).*0;
matWVLST = matWVLST + tmp2.*valDELTIME;

%% Recreating wake point matrix, and regenerating wake HDVE parameters
matWAKEGEOM(:,:,1) = matWVLST(matWDVE(:,1,1),:);
matWAKEGEOM(:,:,2) = matWVLST(matWDVE(:,2,1),:);
matWAKEGEOM(:,:,3) = matWVLST(matWDVE(:,3,1),:);

matWETA = nan(valWNELE, 1);
matWETA(reshape([1:valWSIZE.*2:valWNELE]' + [1:valWSIZE]-1, [], 1),1) = 1;
matWETA(reshape([1:valWSIZE.*2:valWNELE]' + [1:valWSIZE] + valWSIZE - 1, [], 1),1) = 2;

[~, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, ...
    matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG] = fcnTRIANG(matWAKEGEOM, 'WAKE', matWETA);

end

