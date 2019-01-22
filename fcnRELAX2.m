function [matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG, matWAKEGEOM] = ...
                fcnRELAX2(valTIMESTEP, matWELOC, matWEATT, matWEIDX, matUINF, valDELTIME, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, ...
                matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST)
% Yeh's Method       

wdves = [flipud([(1:valWSIZE)+(valWSIZE*2.*(1:valTIMESTEP - 2)')]) flipud(valWSIZE*2.*(2:valTIMESTEP - 1)')];
loc_edge = ones(size(wdves));
loc_edge(:,end) = loc_edge(:,end) + 1;
glob_edge = reshape(matWEIDX(sub2ind(size(matWEIDX),wdves(:),loc_edge(:))),size(wdves)); 

fpg = matWCENTER(wdves,:);
q_ind = fcnSDVEVEL(fpg, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER) + fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER).*valDELTIME;
dir_ind = q_ind./(sqrt(sum(q_ind.^2,2)));

vec_edge = matWVLST(matWELST(glob_edge,2),:) - matWVLST(matWELST(glob_edge,1),:);
dir_edge = vec_edge./(sqrt(sum(vec_edge.^2,2)));

matWVLST(matWELST(glob_edge,2),:) = vec_edge + q_ind;

hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3));
hold off

%% Recreating wake point matrix, and regenerating wake HDVE parameters
matWAKEGEOM(:,:,1) = matWVLST(matWDVE(:,1,1),:);
matWAKEGEOM(:,:,2) = matWVLST(matWDVE(:,2,1),:);
matWAKEGEOM(:,:,3) = matWVLST(matWDVE(:,3,1),:);

[~, matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, ...
    matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG] = fcnTRIANG(matWAKEGEOM);

