function [matWELST, matWVLST, matWDVE, valWNELE, matWEIDX, matWPLEX, matWDVECT, matWVATT, matWCENTER, matWROTANG, matWAKEGEOM, matWVGRID] = ...
    fcnRELAX2(valTIMESTEP, matUINF, valDELTIME, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, ...
    matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST, matWVATT, matWEIDX, vecDVESYM, vecWDVESYM, vecWSYM, matWVGRID, valPRESTEPS, matWPLANE)

% idx = flipud(repmat([1:valWSIZE, valWSIZE*2], valTIMESTEP, 1) + [0:(valWSIZE*2):(valWSIZE*valTIMESTEP*2 - 1)]');
idx = flipud(repmat([(1:valWSIZE) + valWSIZE, valWSIZE*2], valTIMESTEP, 1) + [0:(valWSIZE*2):(valWSIZE*valTIMESTEP*2 - 1)]');
q_ind = nan(valWNELE, 3);
q_ind(idx,:) = fcnSDVEVEL(matWCENTER(idx,:), valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM) + fcnSDVEVEL(matWCENTER(idx,:), valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, vecWDVESYM);

% Project onto plane
for i = 1:size(idx,1)
    q_ind(idx(i,:),:) = q_ind(idx(i,:),:) - dot(q_ind(idx(i,:),:), repmat(matWPLANE(i,:), size(q_ind(idx(i,:),:), 1), 1), 2);
end
idx_for = [idx(1,:); idx(1:end-1,:)];
idx_aft = [idx(2:end,:); idx(end,:)];
q_ind(idx,:) = (q_ind(idx_for,:) + 2.*q_ind(idx,:) + q_ind(idx_aft,:))./4;

% hold on
% quiver3(matWCENTER(idx,1), matWCENTER(idx,2), matWCENTER(idx,3), q_ind(idx,1), q_ind(idx,2), q_ind(idx,3))
% hold off

tmp2 = zeros(size(matWVLST));
for j = 1:valWSIZE + 1
    if j <= valWSIZE
        tmp2(matWDVE(idx(:,j),1),:) = q_ind(idx(:,j),:);
    else
        tmp2(matWDVE(idx(:,j),3),:) = q_ind(idx(:,j),:);
    end
end





%% Moving
% Not moving some vertices
tmp2(matWELST(vecWLE,:),:) = 0; % Trailing edge of wing
tmp2(matWELST(vecWSYM,:),2) = 0; % No y-component on symmetry line
oldest_edge = matWEIDX((1:valWSIZE) + valWSIZE,3); % Trailing edge of oldest wake row
tmp2(matWELST(oldest_edge,:),:) = 0;

% hold on
% quiver3(matWVLST(:,1), matWVLST(:,2), matWVLST(:,3), tmp2(:,1), tmp2(:,2), tmp2(:,3))
% hold off

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


