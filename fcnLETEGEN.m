function [vecLE, vecLEDVE, vecTE, vecTEDVE, vecSYM, vecSYMDVE, matELST, matDVEGRID, vecCHORD, vecSPAN] = fcnLETEGEN(strATYPE, matVLST, matELST, matEATT, matLEPOINTS, matTEPOINTS, vecSYM_old, matEIDX, vecM)

vecLE = [];
vecTE = [];
vecLEDVE = [];
vecTEDVE = [];
vecSYM = [];
vecSYMDVE = [];

if ~isempty(matLEPOINTS)
    [~, idxle(:,1)] = ismember(matLEPOINTS(:,:,1),matVLST,'rows');
    [~, idxle(:,2)] = ismember(matLEPOINTS(:,:,2),matVLST,'rows');
    [~, vecLE] = ismember(sort(idxle,2), sort(matELST,2),'rows');
    vecLEDVE = nonzeros(sort(matEATT(vecLE,:),2,'descend'));
end

if ~isempty(matTEPOINTS) && strcmpi(strATYPE{2},'PANEL')
    [~, idxte(:,1)] = ismember(matTEPOINTS(:,:,1),matVLST,'rows');
    [~, idxte(:,2)] = ismember(matTEPOINTS(:,:,2),matVLST,'rows');
    [~, vecTE] = ismember(sort(idxte,2), sort(matELST,2),'rows');
    vecTEDVE = sort(matEATT(vecTE,:),2,'descend');
elseif ~isempty(matTEPOINTS) && strcmpi(strATYPE{2},'THIN')
    [~, idxte(:,1)] = ismember(matTEPOINTS(:,:,1),matVLST,'rows');
    [~, idxte(:,2)] = ismember(matTEPOINTS(:,:,2),matVLST,'rows');
    [~, vecTE] = ismember(sort(idxte,2), sort(matELST,2), 'rows');
    matELST(vecTE,:) = idxte;
    vecTEDVE = nonzeros(sort(matEATT(vecTE,:),2,'descend'));
end

if any(vecSYM_old)
    symverts = find(matVLST(:,2) == 0);
    vecSYM = find(all(ismember(matELST, symverts),2));
    vecSYMDVE = nonzeros(sort(matEATT(vecSYM,:),2,'descend'));
end

matDVEGRID(1,:) = vecLEDVE;
for i = 2:vecM*2
    [~,matDVEGRID(i,:)] = ismember(matEIDX(matDVEGRID(i-1,:), 3), matEIDX(:,2));
end

le_pt = (matVLST(matELST(vecLE,1),:) + matVLST(matELST(vecLE,2),:))./2;
te_pt = (matVLST(matELST(vecTE,1),:) + matVLST(matELST(vecTE,2),:))./2;

vecCHORD = sqrt(sum((te_pt - le_pt).^2,2));

vecSPAN = sqrt(sum([matVLST(matELST(vecLE,1),:) - matVLST(matELST(vecLE,2),:)].^2,2));

end

