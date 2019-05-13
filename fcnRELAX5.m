function [matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG, matWAKEGEOM] = ...
                fcnRELAX5(valTIMESTEP, matUINF, valDELTIME, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, ...
                matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST, matWVATT, matWEIDX, vecDVESYM, vecWDVESYM, vecWSYM)

dvegrid = flipud(repmat([1:valWSIZE, valWSIZE*2], valTIMESTEP, 1) + [0:(valWSIZE*2):(valWSIZE*valTIMESTEP*2 - 1)]');
edges = reshape([matWEIDX(dvegrid(:,1:end-1), 1)' matWEIDX(dvegrid(:,end), 2)'],size(dvegrid));  
verts(:,:,1) = [reshape(matWELST(edges,1), size(edges))];
verts(:,:,2) = [reshape(matWELST(edges,2), size(edges))];

fpg = (matWVLST(verts(:,:,1),:) + matWVLST(verts(:,:,2),:))./2;
q_ind = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM) + fcnSDVEVEL(fpg, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, vecWDVESYM);

idx = reshape(1:size(fpg,1), size(edges));

% hold on
% % scatter3(fpg(:,1), fpg(:,2), fpg(:,3), 'o');
% quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3));
% % scatter3(matWCENTER(dvegrid(:),1),matWCENTER(dvegrid(:),2),matWCENTER(dvegrid(:),3), '^')
% hold off

%% Moving wake vertices

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

tmp2(matWELST(vecWLE,:),:) = tmp2(matWELST(vecWLE,:),:).*0; % Trailing edge of wing
tmp2(matWELST(vecWTE,:),:) = tmp2(matWELST(vecWTE,:),:).*0; % Trailing edge of LE wake element row
oldest_edge = matWEIDX((1:valWSIZE) + valWSIZE,3); % Trailing edge of oldest wake row
tmp2(matWELST(vecWSYM,:),2) = tmp2(matWELST(vecWSYM,:),2).*0; % No y-component on symmetry line
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

