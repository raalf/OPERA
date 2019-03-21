function [matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG, matWAKEGEOM] = ...
                fcnRELAX(valTIMESTEP, matUINF, valDELTIME, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, ...
                matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST, matWVATT, matWEIDX, vecDVESYM, vecWDVESYM, vecWSYM)

% Finding induced velocities at all wake vertices
q_ind = fcnSDVEVEL(matWCENTER, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM) + fcnSDVEVEL(matWCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, vecWDVESYM);

%% Moving wake vertices
matWCENTER = matWCENTER + q_ind.*valDELTIME;

oldest_edge = matWEIDX((1:valWSIZE) + valWSIZE,3);
% Not moving trailing edge of wing, trailing edge of newest wake row,
% trailing edge of oldest wake row
dont_move = unique([matWELST(oldest_edge,:); matWELST(vecWLE,:); matWELST(vecWTE,:)]);

move_em = setdiff([1:size(matWVLST,1)]', dont_move);
if any(move_em)
    for i = 1:size(move_em,1)
        vtx = move_em(i);
        
        dves = matWVATT(vtx, ~isnan(matWVATT(vtx,:)));
        weights = sqrt(sum((matWCENTER(dves',:) - matWVLST(vtx,:)).^2, 2)); 
        weights = weights./sum(weights);
        matWVLST(vtx,:) = matWVLST(vtx,:) + sum(repmat(weights, 1, 3).*q_ind(dves',:),1).*valDELTIME;    
    end
end


% hFig1 = gcf;
% fcnPLOTWAKE(1, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, valWSIZE, 1)

% for i = 1:size(matWVATT,1)
%    tmp2(i,:) = mean(q_ind(~isnan(matWVATT(i,:)),:),1);
% end

% tmp2(matWELST(vecWLE,:),:) = tmp2(matWELST(vecWLE,:),:).*0; % Trailing edge of wing
% tmp2(matWELST(vecWTE,:),:) = tmp2(matWELST(vecWTE,:),:).*0; % Trailing edge of LE wake element row
% oldest_edge = matWEIDX((1:valWSIZE) + valWSIZE,3); % Trailing edge of oldest wake row
% tmp2(matWELST(oldest_edge,:),:) = tmp2(matWELST(oldest_edge,:),:).*0;
% matWVLST = matWVLST + tmp2.*valDELTIME;

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

