function [TR, ADJE, ELST, VLST, DVE, NELE, EATT, EIDX, ELOC, DNORM, PLEX, DVECT, ALIGN, VATT] = fcnTRIANG(STL, ATYPE)
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
%   DNORM - Matrix of face normals for each DVE. NELE x 3 matrix (columns are (x,y,z))
%   PLEX - Matrix of local eta-xi coordinates of vertices 
%   DVECT - Normal vectors
%   ALIGN - Gives alignment of eta and xi axis between adjacent HDVEs in EATT rows. ALIGN(:,:,1) is alignment of eta in terms of eta and xi of the adjacent HDVE
            % ALIGN(:,:,2) is alignment of xi in terms of eta and xi in adjacent HDVE.
% This function was written to work with the panel code, and lifting surface with and without a split (of SDEG 3). Modelling
% a split in the wing will result in the ADJE matrix having a third dimension. It could likely easily be modified to accept
% a split of SDEG > 3.
% T.D.K 2016-08-15. 212-230 KING ST E, TORONTO, ONTARIO, CANADA, M5A-1K5

[temp, ~] = fcnSTLREAD(STL);

% Number of DVEs
NELE = size(temp(:,:,1),1);

% Getting unique vertices, and switching from (x,y,z) to indices with reference to master list VLST
[VLST,~,j] = unique([temp(:,:,1); temp(:,:,2); temp(:,:,3)],'rows','stable');
DVE(:,:,1) = reshape(j,NELE,3);

if strcmp(ATYPE,'LS')
    % Removing duplicate triangles (as the STL is 2 surfaces ontop of each other)
    DVE(:,:,1) = sort(DVE(:,:,1),2);
    [~, ia, ~] = unique(DVE(:,:,1),'rows');
    DVE = DVE(ia,:,:);
    NELE = NELE/2;
end

% Converting above data to triangulation
TR = triangulation(DVE(:,:,1),VLST);

ELST = edges(TR); % List of unique edges
DNORM = -faceNormal(TR);

DVE(:,:,2) = faceNormal(TR); % Normal
DVE(:,:,3) = incenter(TR); % incenter of triangle
DVE(:,:,4) = circumcenter(TR);

VATT = vertexAttachments(TR);

%% Finding edge attachement matrix (which DVEs share which edge)

% Different numbers of edge attachements means we can't do a straight cell2mat on edgeAttachements(TR,ELST) if using VAP2
temp2 = edgeAttachments(TR, ELST);

cellsz = cell2mat(cellfun(@size,temp2,'uni',false));
idx = cellsz(:,2);
if max(idx) == 3
    temp3(idx==3,:) = cell2mat(temp2(idx==3,:));
    temp3(idx==2,2:3) = cell2mat(temp2(idx==2,:));
else
    temp3(idx==2,1:2) = cell2mat(temp2(idx==2,:));
end
temp3(idx==1,1) = cell2mat(temp2(idx==1,:));

EATT = sort(temp3,2); % List of unique edge attachements (by element #), and sort it on each row

%% Mapping global edge number to local edge number in EIDX
% This may be improved to fewer lines

% Making an NELE*3 x 2 matrix of vertices between which are edges. The first third is for the first edge of each HDVE,
% the second third is for the second edge of each HDVE, etc.
dve1 = reshape([DVE(:,1:2,1) DVE(:,2:3,1) DVE(:,3,1) DVE(:,1,1)],NELE,2,3);
dve2 = [dve1(:,:,1); dve1(:,:,2); dve1(:,:,3)];

% Finding the number of edges
nedg = length(ELST(:,1));

% Finding which edge occurs on which DVE, returned as a cell array
% The reverse is searched for as well to find both instances of an edge
[~,a1] = ismembertol(ELST,dve2,'ByRows',true,'OutputAllIndices',true);
[~,b1] = ismembertol([ELST(:,2) ELST(:,1)],dve2,'ByRows',true,'OutputAllIndices',true);

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
EIDX = temp4(idx2);
% As it was an nedg*2 x 1 list, we remove nedg for values > nedg to start at 1 again
temp8 = (floor(EIDX./nedg).*nedg);
EIDX(EIDX>nedg) = EIDX(EIDX>nedg) - temp8(EIDX>nedg);
EIDX(EIDX==0) = nedg;
EIDX = reshape(EIDX,[],3);

%% Computing ELOC

[i3,~,~] = find(EATT>0); % Finding indices of nonzero elements in EATT
[i,j,~] = find(EIDX(EATT(EATT>0),:) == repmat(i3,1,3)); % Finding where these indices are equal to whats in EIDX

temp31 = sortrows([i j]); % Sorting by column 1 so column 2 is the correct order
ELOC(EATT>0,1) = temp31(:,2); % ELOC is column 2, indexing it appropriately
ELOC = reshape(ELOC,nedg,[]);

%% Computing the adjaceny matrix (which DVEs are adjacent to which, columns are sensitive to edge number)

SDEG = length(EATT(1,:)); % The degree of the split

% Neighbouring elements along the first local edge of each HDVE
temp10 = EATT(EIDX(:,1),:); % Getting all HDVEs that share the local first edge of every HDVE
temp10(temp10==repmat([1:NELE]',1,SDEG)) = 0; % Removing HDVE from its own adjacency
temp10 = sort(temp10,2); % Sorting to sift out the zeros from the columns
temp10(:,1) = []; % Removing extra column of zeros after the sort
ADJE(:,1,:) = reshape(temp10,NELE,1,SDEG-1); % Placing values in ADJE matrix

clear temp10

% Local edge #2
temp10 = EATT(EIDX(:,2),:); % Getting all HDVEs that share the local first edge of every HDVE
temp10(temp10==repmat([1:NELE]',1,SDEG)) = 0; % Removing HDVE from its own adjacency
temp10 = sort(temp10,2); % Sorting to sift out the zeros from the columns
temp10(:,1) = []; % Removing extra column of zeros after the sort
ADJE(:,2,:) = reshape(temp10,NELE,1,SDEG-1); % Placing values in ADJE matrix

clear temp10

% Local edge #3
temp10 = EATT(EIDX(:,3),:); % Getting all HDVEs that share the local first edge of every HDVE
temp10(temp10==repmat([1:NELE]',1,SDEG)) = 0; % Removing HDVE from its own adjacency
temp10 = sort(temp10,2); % Sorting to sift out the zeros from the columns
temp10(:,1) = []; % Removing extra column of zeros after the sort
ADJE(:,3,:) = reshape(temp10,NELE,1,SDEG-1); % Placing values in ADJE matrix

clear temp10

% If missing an adjacency with panel code, then there is a discontinuity in the geometry somewhere
if isempty(find(isnan(ADJE))) == 0 && strcmp(ATYPE,'PC')
    disp('Problem with geometry in fcnTRIANG.')
end

%% Local HDVE Eta-Xi Axis

P = permute(reshape(TR.Points(TR.ConnectivityList',:)',3,3,[]),[2 1 3]);
[PLEX, DVECT] = fcnTRITOLEX(P, DNORM);

% % Plotting global and local to visualize
% test_num = 396;
% hFig2 = figure(2);
% clf(2);
% subplot(2,1,1)
% patch(VLEX(:,1,test_num),VLEX(:,2,test_num),'b')
% alpha(0.5);
% xlabel('u-direction','FontSize',15);
% ylabel('v-direction','FontSize',15);
% axis equal
% grid on
% box on
% subplot(2,1,2)
% patch(P(:,1,test_num),P(:,2,test_num),P(:,3,test_num),'r')
% alpha(0.5);
% xlabel('X-direction','FontSize',15);
% ylabel('Y-direction','FontSize',15);
% zlabel('Z-direction','FontSize',15);
% axis equal
% grid on
% box on

%% Finding DVE Alignment Matrix ALIGN

ALIGN = repmat(zeros(size(EATT)),1,1,2);

% Projecting vector from DVE onto adjacent DVE
idx = all(EATT,2); % All edges that split 2 DVEs

% Eta direction of DVE projected onto adjacent DVE
vec1 = DVECT(EATT(idx,1),:,1) - repmat((dot(DVECT(EATT(idx,1),:,1),DVECT(EATT(idx,2),:,3), 2)./(sqrt(sum(abs(DVECT(EATT(idx,2),:,3).^2),2)).^2)),1,3).*DVECT(EATT(idx,2),:,3);

% Xi direction of DVE projected onto adjacent DVE
vec2 = DVECT(EATT(idx,1),:,2) - repmat((dot(DVECT(EATT(idx,1),:,2),DVECT(EATT(idx,2),:,3), 2)./(sqrt(sum(abs(DVECT(EATT(idx,2),:,3).^2),2)).^2)),1,3).*DVECT(EATT(idx,2),:,3);

% If any of the projected vectors are equal to the normal, then take the cross product of the other projected vector and the normal
idx1 = ~any(vec1,2);
vec1(idx1,:) = cross(DVECT(EATT(idx1,2),:,3),vec2(idx1,:)); 
idx2 = ~any(vec2,2);
vec2(idx2,:) = cross(vec1(idx2,:),DVECT(EATT(idx2,2),:,3)); 

% Finding local eta, xi of projections on adjacent DVE
ALIGN(idx,1,1) = dot(vec1,DVECT(EATT(idx,2),:,1),2);
ALIGN(idx,2,1) = dot(vec1,DVECT(EATT(idx,2),:,2),2);

ALIGN(idx,1,2) = dot(vec2,DVECT(EATT(idx,2),:,1),2);
ALIGN(idx,2,2) = dot(vec2,DVECT(EATT(idx,2),:,2),2);


clearvars -except TR ADJE ELST VLST DVE NELE EATT EIDX ELOC DNORM PLEX DVECT ALIGN VATT

end

