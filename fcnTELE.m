function [vecTE, vecLE] = fcnTELE(matTEPOINTS, matLEPOINTS, matVLST, matELST)
% Identifies which edge (by edge ID) is a trailing or leading edge

% Trailing Edge
[~, vertex1] = ismember(matTEPOINTS(:,:,1), matVLST, 'rows');
[~, vertex2] = ismember(matTEPOINTS(:,:,2), matVLST, 'rows');

[~, vecTE] = ismember(sort([vertex1 vertex2],2), sort(matELST,2),'rows');

% Leading Edge
[~, vertex1] = ismember(matLEPOINTS(:,:,1), matVLST, 'rows');
[~, vertex2] = ismember(matLEPOINTS(:,:,2), matVLST, 'rows');

[~, vecLE] = ismember(sort([vertex1 vertex2],2), sort(matELST,2),'rows');

end

