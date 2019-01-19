function [vecLE, vecLEDVE, vecTE, vecTEDVE, matSPANDIR] = fcnLETEGEN(strATYPE, valNELE, matVLST, matELST, matDVECT, matEATT, matLEPOINTS, matTEPOINTS)

matSPANDIR = repmat([0 1 0],valNELE,1);
matSPANDIR = matSPANDIR - (dot(matSPANDIR, matDVECT(:,:,3),2)).*matDVECT(:,:,3);
vecLE = [];
vecTE = [];
vecLEDVE = [];
vecTEDVE = [];

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
    vecTEDVE = nonzeros(sort(matEATT(vecTE,:),2,'descend'));
end

end

