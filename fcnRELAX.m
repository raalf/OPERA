function [matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG, matWAKEGEOM] = ...
                fcnRELAX(vecUINF, valDELTIME, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, ...
                matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST)
            
% Relaxes wake by moving points in the wake vertex list matWVLST, and updating the local vectors in matWDVECT

%% Finding induced velocities at all wake vertices

% for i = 1:size(matWVLST,1)
%     [attached_dves,~] = find(matWDVE == i);
%     fpg = matWCENTER(attached_dves,:);
%     q_ind(i,:) = mean(fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER) + fcnSDVEVEL(fpg, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER), 1);
% end

q_ind = fcnSDVEVEL(matWVLST+[0.1 0 0], valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER) + fcnSDVEVEL(matWVLST, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER);

hold on
quiver3(matWVLST(:,1), matWVLST(:,2), matWVLST(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3),'r')
hold off

%% Identifying trailing edges, we won't relax these points

% Finding vertices common between wing and wake, these are the trailing edge vertices
% dontmove = [matWVLST(matWELST(vecWLE,1),:); matWVLST(matWELST(vecWLE(end),2),:); matWVLST(matWELST(vecWTE,1),:); matWVLST(matWELST(vecWTE(end),2),:)];
oldest_dve = [1:(valWSIZE*2)]';
dontmove = [matWVLST(matWELST(vecWLE,1),:); matWVLST(matWELST(vecWLE(end),2),:); matWVLST(matWDVE(oldest_dve,:),:)];
[idxTEV, ~] = ismember(matWVLST, dontmove, 'rows');

% Setting induced velocities at these points to zero
q_ind(idxTEV,:) = repmat([0 0 0], length(nonzeros(idxTEV)), 1);

%% Moving wake vertices
matWVLST = matWVLST + q_ind.*valDELTIME;

%% Recreating wake point matrix, and regenerating wake HDVE parameters
matWAKEGEOM(:,:,1) = matWVLST(matWDVE(:,1,1),:);
matWAKEGEOM(:,:,2) = matWVLST(matWDVE(:,2,1),:);
matWAKEGEOM(:,:,3) = matWVLST(matWDVE(:,3,1),:);

[~, matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, ...
    matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG] = fcnTRIANG(matWAKEGEOM);

