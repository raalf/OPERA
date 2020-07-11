function [vecLE, vecLEDVE, vecTE, vecTEDVE, matELST, vecCHORD, vecSPAN] = fcnLETEGEN(matVLST, matELST, matEATT, matLEPOINTS, matTEPOINTS)

vecLE = [];
vecLEDVE = [];

if ~isempty(matLEPOINTS)
    [~, idxle(:,1)] = ismember(matLEPOINTS(:,:,1),matVLST,'rows');
    [~, idxle(:,2)] = ismember(matLEPOINTS(:,:,2),matVLST,'rows');
    [~, vecLE] = ismember(sort(idxle,2), sort(matELST,2),'rows');
    vecLEDVE = nonzeros(sort(matEATT(vecLE,:),2,'descend'));
end

[~, idxte(:,1)] = ismember(matTEPOINTS(:,:,1),matVLST,'rows');
[~, idxte(:,2)] = ismember(matTEPOINTS(:,:,2),matVLST,'rows');
[~, vecTE] = ismember(sort(idxte,2), sort(matELST,2), 'rows');
matELST(vecTE,:) = idxte;
vecTEDVE = nonzeros(sort(matEATT(vecTE,:),2,'descend'));

le_pt = (matVLST(matELST(vecLE,1),:) + matVLST(matELST(vecLE,2),:))./2;
te_pt = (matVLST(matELST(vecTE,1),:) + matVLST(matELST(vecTE,2),:))./2;

vecCHORD = sqrt(sum((te_pt - le_pt).^2,2));

vecSPAN = sqrt(sum([matVLST(matELST(vecLE,1),:) - matVLST(matELST(vecLE,2),:)].^2,2));

end

