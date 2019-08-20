function [vecLE, vecLEDVE, vecTE, vecTEDVE, vecSYM, vecSYMDVE, matELST] = fcnLETEGEN(strATYPE, matVLST, matELST, matEATT, matLEPOINTS, matTEPOINTS, vecSYM_old)

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



end

