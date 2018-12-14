function [TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
            matELOC, matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG, matCONTROL] = fcnTRIANG(POINTS)
% This function reads the STL and creates the HDVE matrices.
% Inputs:
%   POINTS - n x 3 x 3 matrix, where columns are (x,y,z) and depth is vertex number
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
%   VATT - NELE x ? matrix of which elements are attached to which vertices
%   VNORM - NELE x 3 matrix of the averaged normals of all elements attached to a vertex (used for flow tangency)
% This function was written to work with the panel code, and lifting surface with and without a split (of SDEG 3). Modelling
% a split in the wing will result in the ADJE matrix having a third dimension. It could likely easily be modified to accept
% a split of SDEG > 3.
% T.D.K 2016-08-15. 212-230 KING ST E, TORONTO, ONTARIO, CANADA, M5A-1K5

% Number of DVEs
valNELE = size(POINTS(:,:,1),1);

% Getting unique vertices, and switching from (x,y,z) to indices with reference to master list VLST
% [matVLST,~,j] = unique([POINTS(:,:,1); POINTS(:,:,2); POINTS(:,:,3)],'rows','stable');
tol = 1e-14;
[matVLST,~,j] = uniquetol([POINTS(:,:,1); POINTS(:,:,2); POINTS(:,:,3)], tol, 'ByRows',true);
matDVE(:,:,1) = reshape(j,[],3);

% % Removing duplicate elements from OpenVSP STL when infinitely thin
% if strcmp(strATYPE,'LS')
%    [~, ia, ~] = unique(sort(matDVE(:,:,1),2),'rows','stable');
%    matDVE = matDVE(ia,:,:);
%    valNELE = length(matDVE(:,1,1)); 
% end

% Converting above data to triangulation
TR = triangulation(matDVE(:,:,1),matVLST);

matELST = edges(TR); % List of unique edges
DNORM = -faceNormal(TR);
% DNORM(1,:) = [0 0 1] %TEMPORARY DON'T KEEP

% matDVE(:,:,2) = faceNormal(TR); % Normal
% matCENTER = incenter(TR); % incenter of triangle
matCENTER = (matVLST(matDVE(:,1),:) + matVLST(matDVE(:,2),:) + matVLST(matDVE(:,3),:))./3;
matCONTROL = matCENTER;
% matCONTROL = (matVLST(matDVE(:,2),:) + matVLST(matDVE(:,1),:) + matCENTER)./3;

%% Finding edge attachement matrix (which DVEs share which edge)

temp2 = edgeAttachments(TR, matELST);

% Different numbers of edge attachements means we can't do a straight cell2mat on edgeAttachements(TR,ELST)
cellsz = cell2mat(cellfun(@size,temp2,'uni',false));
idx = cellsz(:,2);
if max(idx) == 3
    temp3(idx==3,:) = cell2mat(temp2(idx==3,:));
    temp3(idx==2,2:3) = cell2mat(temp2(idx==2,:));
else
    temp3(idx==2,1:2) = cell2mat(temp2(idx==2,:));
end
temp3(idx==1,1) = cell2mat(temp2(idx==1,:));

matEATT = sort(temp3,2); % List of unique edge attachements (by element #), and sort it on each row

%% Mapping global edge number to local edge number in EIDX
% This may be improved to fewer lines

% Making an NELE*3 x 2 matrix of vertices between which are edges. The first third is for the first edge of each HDVE,
% the second third is for the second edge of each HDVE, etc.
dve1 = reshape([matDVE(:,1:2,1) matDVE(:,2:3,1) matDVE(:,3,1) matDVE(:,1,1)],valNELE,2,3);
dve2 = [dve1(:,:,1); dve1(:,:,2); dve1(:,:,3)];

% Finding the number of edges
nedg = length(matELST(:,1));

% Finding which edge occurs on which DVE, returned as a cell array
% The reverse is searched for as well to find both instances of an edge
[~,a1] = ismembertol(matELST,dve2,'ByRows',true,'OutputAllIndices',true);
[~,b1] = ismembertol([matELST(:,2) matELST(:,1)],dve2,'ByRows',true,'OutputAllIndices',true);

clear temp2 temp3

% Reorganizing edge numbers into a list of which edge occurs on which DVE (an edge can appear maximum of max(idx))
cellsza = cell2mat(cellfun(@size,a1,'uni',false));
idxa = cellsza(:,1);
cellszb = cell2mat(cellfun(@size,b1,'uni',false));
idxb = cellszb(:,1);

if max(idxa) == 3 || max(idxb) == 3
    temp5(idxa==3,:) = reshape(cell2mat(a1(idxa==3,:)),3,[])';
    temp6(idxb==3,:) = reshape(cell2mat(b1(idxb==3,:)),3,[])';
end
temp5(idxa==2,1:2) = reshape(cell2mat(a1(idxa==2,:)),2,[])';
temp6(idxb==2,1:2) = reshape(cell2mat(b1(idxb==2,:)),2,[])';
temp5(idxa==1,1) = cell2mat(a1(idxa==1,:));
temp6(idxb==1,1) = cell2mat(b1(idxb==1,:));
temp7 = [temp5 temp6];

% Putting them as a nedg*2 x 1 list
temp4 = find(reshape(temp7,[],1)>0);
% Finding the indices that order the above list appropriately
[~,idx2] = sort(temp7(find(temp7>0)));
matEIDX = temp4(idx2);
% As it was an nedg*2 x 1 list, we remove nedg for values > nedg to start at 1 again
temp8 = (floor(matEIDX./nedg).*nedg);
matEIDX(matEIDX>nedg) = matEIDX(matEIDX>nedg) - temp8(matEIDX>nedg);
matEIDX(matEIDX==0) = nedg;
matEIDX = reshape(matEIDX,[],3);

%% Computing ELOC
[i3,~,~] = find(matEATT>0); % Finding indices of nonzero elements in EATT
[i,j,~] = find(matEIDX(matEATT(matEATT>0),:) == repmat(i3,1,3)); % Finding where these indices are equal to whats in EIDX

temp31 = sortrows([i j]); % Sorting by column 1 so column 2 is the correct order
matELOC(matEATT>0,1) = temp31(:,2); % ELOC is column 2, indexing it appropriately
matELOC = reshape(matELOC,nedg,[]);

%% Computing the adjaceny matrix (which DVEs are adjacent to which, columns are sensitive to edge number)
SDEG = length(matEATT(1,:)); % The degree of the split

% Neighbouring elements along the first local edge of each HDVE
temp10 = matEATT(matEIDX(:,1),:); % Getting all HDVEs that share the local first edge of every HDVE
temp10(temp10==repmat([1:valNELE]',1,SDEG)) = 0; % Removing HDVE from its own adjacency
temp10 = sort(temp10,2); % Sorting to sift out the zeros from the columns
temp10(:,1) = []; % Removing extra column of zeros after the sort
matADJE(:,1,:) = reshape(temp10,valNELE,1,SDEG-1); % Placing values in ADJE matrix

clear temp10

% Local edge #2
temp10 = matEATT(matEIDX(:,2),:); % Getting all HDVEs that share the local first edge of every HDVE
temp10(temp10==repmat([1:valNELE]',1,SDEG)) = 0; % Removing HDVE from its own adjacency
temp10 = sort(temp10,2); % Sorting to sift out the zeros from the columns
temp10(:,1) = []; % Removing extra column of zeros after the sort
matADJE(:,2,:) = reshape(temp10,valNELE,1,SDEG-1); % Placing values in ADJE matrix

clear temp10

% Local edge #3
temp10 = matEATT(matEIDX(:,3),:); % Getting all HDVEs that share the local first edge of every HDVE
temp10(temp10==repmat([1:valNELE]',1,SDEG)) = 0; % Removing HDVE from its own adjacency
temp10 = sort(temp10,2); % Sorting to sift out the zeros from the columns
temp10(:,1) = []; % Removing extra column of zeros after the sort
matADJE(:,3,:) = reshape(temp10,valNELE,1,SDEG-1); % Placing values in ADJE matrix

clear temp10

% If missing an adjacency with panel code, then there is a discontinuity in the geometry somewhere
if isempty(find(isnan(matADJE))) == 0 && strcmp(ATYPE,'PC')
    disp('Problem with geometry in fcnTRIANG.')
end

%% Local HDVE Xi-eta Axis
P = permute(reshape(TR.Points(TR.ConnectivityList',:)',3,3,[]),[2 1 3]);
[matPLEX, matDVECT, matROTANG] = fcnTRITOLEX(P, DNORM, matCENTER);


%% Vertex attachements and normal averages
matVATT = vertexAttachments(TR);

% Turning the cell array VATT into a matrix - A.Y. Method 2016-09-20
cellsza = cell2mat(cellfun(@size,matVATT,'uni',false));
idxa = cellsza(:,2);
idx2 = num2cell(max(cell2mat(cellfun(@length,matVATT,'uni',false)))-idxa);
matVATT = cell2mat(cellfun(@(x,y) padarray(x,[0 y],NaN,'post'), matVATT, idx2, 'uni', false));

% Finding the averaged face-normals of all elements attached to the vertices
idx50 = matVATT>0;
temp50 = permute(matDVECT(matVATT(idx50),:,3),[1,3,2]);
temp51 = zeros(size(matVATT));
temp52 = zeros(size(matVATT));
temp53 = zeros(size(matVATT));
temp51(idx50) = temp50(:,:,1);
temp52(idx50) = temp50(:,:,2);
temp53(idx50) = temp50(:,:,3);

% Averaging the normals of the x, y, z components of the normals attached to each vertex
matVNORM = [mean(temp51,2) mean(temp52,2) mean(temp53,2)];

% Normalizing these vectors
matVNORM = matVNORM./repmat(sqrt(sum(abs(matVNORM).^2,2)), 1,3);

% clearvars -except TR ADJE ELST VLST DVE NELE EATT EIDX ELOC PLEX DVECT ALIGN VATT VNORM CENTER

end

