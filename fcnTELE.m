function [vecTE, vecLE] = fcnTELE(matTEPOINTS, matLEPOINTS, matVLST, matELST)
% Identifies which edge (by edge ID) is a trailing or leading edge

% Trailing Edge
% [~, vertex1] = ismember(matTEPOINTS(:,:,1), matVLST, 'rows');
% [~, vertex2] = ismember(matTEPOINTS(:,:,2), matVLST, 'rows');
% [~, vecTE] = ismember(sort([vertex1 vertex2],2), sort(matELST,2),'rows');

tol = 1e-6;

[~, vertex1] = ismembertol(matTEPOINTS(:,:,1), matVLST, tol, 'ByRows', true);
[~, vertex2] = ismembertol(matTEPOINTS(:,:,2), matVLST, tol, 'ByRows', true);
[~, vecTE] = ismembertol(sort([vertex1 vertex2],2), sort(matELST,2), tol, 'ByRows', true);

% Leading Edge
% [~, vertex1] = ismember(matLEPOINTS(:,:,1), matVLST, 'rows');
% [~, vertex2] = ismember(matLEPOINTS(:,:,2), matVLST, 'rows');
% [~, vecLE] = ismember(sort([vertex1 vertex2],2), sort(matELST,2),'rows');

[~, vertex1] = ismembertol(matLEPOINTS(:,:,1), matVLST, tol, 'ByRows', true);
[~, vertex2] = ismembertol(matLEPOINTS(:,:,2), matVLST, tol, 'ByRows', true);
[~, vecLE] = ismembertol(sort([vertex1 vertex2],2), sort(matELST,2), tol, 'ByRows', true);

end

