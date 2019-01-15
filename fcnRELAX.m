function [matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG, matWAKEGEOM] = ...
                fcnRELAX(vecUINF, valDELTIME, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, ...
                matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST)
            
% Relaxes wake by moving points in the wake vertex list matWVLST, and updating the local vectors in matWDVECT

%% Finding induced velocities at all wake vertices
s_ind = fcnSDVEVEL(matWVLST, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);
w_ind = fcnSDVEVEL(matWVLST, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER);

hold on
quiver3(matWVLST(:,1), matWVLST(:,2), matWVLST(:,3), w_ind(:,1), w_ind(:,2), w_ind(:,3),'g')
quiver3(matWVLST(:,1), matWVLST(:,2), matWVLST(:,3), s_ind(:,1), s_ind(:,2), s_ind(:,3),'r')
hold off

% s_ind = s_ind.*0;
% w_ind = w_ind.*0;

% for i = 1:valWNELE
%     q_loc = -0.5.*[matWCOEFF(i,3).*matWPLEX(:,1,i) + matWCOEFF(i,4), matWCOEFF(i,1).*matWPLEX(:,2,i) + matWCOEFF(i,2), matWPLEX(:,2,i).*0];
%     
% end


%% Identifying trailing edges, we won't relax these points

% Finding vertices common between wing and wake, these are the trailing edge vertices
% dontmove = [matWVLST(matWELST(vecWLE,1),:); matWVLST(matWELST(vecWLE(end),2),:); matWVLST(matWELST(vecWTE,1),:); matWVLST(matWELST(vecWTE(end),2),:)];
oldest_dve = [1:(valWSIZE*20)]';
dontmove = [matWVLST(matWELST(vecWLE,1),:); matWVLST(matWELST(vecWLE(end),2),:); matWVLST(matWDVE(oldest_dve,:),:)];
[idxTEV, ~] = ismember(matWVLST, dontmove, 'rows');

% Setting induced velocities at these points to zero
s_ind(idxTEV,:) = repmat([0 0 0], length(nonzeros(idxTEV)), 1);
w_ind(idxTEV,:) = repmat([0 0 0], length(nonzeros(idxTEV)), 1);

%% Moving wake vertices
% matWVLST = matWVLST + (s_ind + w_ind + repmat(vecUINF, length(s_ind(:,1)),1)).*valDELTIME; % Needs more work for rotors (assuming freestream magnitude is 1)
matWVLST = matWVLST + (s_ind + w_ind).*valDELTIME;

%% Recreating wake point matrix, and regenerating wake HDVE parameters
matWAKEGEOM(:,:,1) = matWVLST(matWDVE(:,1,1),:);
matWAKEGEOM(:,:,2) = matWVLST(matWDVE(:,2,1),:);
matWAKEGEOM(:,:,3) = matWVLST(matWDVE(:,3,1),:);

[~, matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, ...
    matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG] = fcnTRIANG(matWAKEGEOM);

