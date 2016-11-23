function [TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matELOC, ...
            matPLEX, matDVECT, matALIGN, matVATT, matVNORM, matCENTER] = fcnIMPORTGEOM(strSTL, strATYPE)
% This function reads the STL and creates the HDVE matrices.
% Inputs:
%   STL - .stl filename, string.
%   ATYPE - Analysis type, string. 'LS' - lifting surface, 'PC' - panel code.
% Outputs:
%   TR - MATLAB triangulation
%   ADJE - NELE x 3 x SDEG matrix of adjacent HDVEs, where the row is the HDVE and the columns are the local edge numbers.
%   ELST - List of unique edges, and the vertices they connect.
%   VLST - List of unique vertices, and their global (X,Y,Z) coordinates.
%   DVE - Multi-level matrix. First level is indices from VLST identifying vertices of each HDVE. Second is normal vector.
%   NELE - Number of HDVEs.
%   EATT - Matrix of edge attachements. Rows are unique edges, columns are the HDVEs that use that edge.
%   EIDX - Matrix denoting which edge number is which local edge number. Rows are HDVEs. Columns are local edge numbers 1-3.
%   ELOC - Identifies the local edge number from EIDX in EATT.
%   SDEG - Currently unused. Split degree. No split in the wing has degree 2, as the maximum number of HDVEs per edge is 2.
%   PLEX - Matrix of local eta-xi coordinates of vertices 
%   DVECT - Normal vectors
%   ALIGN - Gives alignment of eta and xi axis between adjacent HDVEs in EATT rows. ALIGN(:,:,1) is alignment of eta in terms of eta and xi of the adjacent HDVE
            % ALIGN(:,:,2) is alignment of xi in terms of eta and xi in adjacent HDVE.
%   VATT - NELE x ? matrix of which elements are attached to which vertices
%   VNORM - NELE x 3 matrix of the averaged normals of all elements attached to a vertex (used for flow tangency)
% This function was written to work with the panel code, and lifting surface with and without a split (of SDEG 3). Modelling
% a split in the wing will result in the ADJE matrix having a third dimension. It could likely easily be modified to accept
% a split of SDEG > 3.
% T.D.K 2016-09-23. KHE 33, 350 VICTORIA STREET, TORONTO, ONTARIO, CANADA M5B-2K3

[temp, ~] = fcnSTLREAD(strSTL);

% Removing duplicate elements from OpenVSP STL when infinitely thin
% if strcmp(strATYPE,'LS')
%     [temp2, idx, ~] = unique(temp,'rows','stable');
% %     temp2(:,:,2:3) = temp(idx,:,2:3);
% %    temp((end/2) + 1:end,:,:) = [];
% end

[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matELOC, matPLEX, matDVECT, matALIGN, matVATT, matVNORM, matCENTER] = fcnTRIANG(strATYPE, temp);

end

