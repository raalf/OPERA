function [matWELST, matWVLST, matWDVE, valWNELE, matWEIDX, matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG, matWAKEGEOM, matWVGRID] = ...
                fcnRELAX2(valTIMESTEP, matUINF, valDELTIME, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, ...
                matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST, matWVATT, matWEIDX, vecDVESYM, vecWDVESYM, vecWSYM, matWVGRID, valPRESTEPS)
            
verts(:,:,1) = matWVGRID(1:end-1,:);
verts(:,:,2) = matWVGRID(2:end,:);

fpg = (matWVLST(verts(:,:,1),:) + matWVLST(verts(:,:,2),:))./2;
q_ind = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM) + fcnSDVEVEL(fpg, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, vecWDVESYM);

idx = reshape(1:size(fpg,1), size(verts(:,:,1)));

edges_for = [idx(1,:); idx(1:end-1,:)];
edges_aft = [idx(2:end,:); idx(end,:)];
q_ind(idx(:),:) = (q_ind(edges_for(:),:) + 2.*q_ind(idx(:),:) + q_ind(edges_aft(:),:))./4;

%% Moving wake vertices
% Getting velocities at all wake vertices
tmp2 = zeros(size(matWVLST));
numvert = zeros(size(matWVLST,1),1);
for i = 1:size(verts,1)
   tv = verts(i,:,:);
   tmp2(tv(:,:,1),:) = tmp2(tv(:,:,1),:) + q_ind(idx(i,:),:);
   tmp2(tv(:,:,2),:) = tmp2(tv(:,:,2),:) + q_ind(idx(i,:),:);
    
   numvert(tv(:,:,1),:) = numvert(tv(:,:,1),:) + 1;
   numvert(tv(:,:,2),:) = numvert(tv(:,:,2),:) + 1;
end
tmp2 = tmp2./numvert;

% Not moving some vertices
tmp2(matWELST(vecWLE,:),:) = 0; % Trailing edge of wing
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

[~, ~, ~, ~, ~, ~, ~, ~, matWPLEX, ...
    matWDVECT, ~, matWVNORM, matWCENTER, matWROTANG] = fcnTRIANG(matWAKEGEOM, 'WAKE', matWETA);

end


