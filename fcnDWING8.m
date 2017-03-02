function [D] = fcnDWING8(strATYPE, matEATT, matPLEX, valNELE, matELOC, matELST, matALIGN, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecSYM, matVATT)

lambda_mid = [ ...
    0.5 0.5 0; ... % Edge 1 mid-point
    0 0.5 0.5; ... % Edge 2 mid-point
    0.5 0 0.5; ... % Edge 3 mid-point
    ];

lambda_vert = [ ...
    1 0 0; ...
    0 1 0; ...
    0 0 1; ...
    ];

% D = sparse(N*5, N*5);

% Making D-matrix a little long, in case we are over-constrained
D = zeros(valNELE*9, valNELE*5);

%% Circulation equations between elements

% Evaluated at the mid-point of each edge which splits two HDVEs
idx = all(matEATT,2); % All edges that split 2 DVEs
nedg = length(matEATT(idx,1));

circ = fcnDCIRC(idx, nedg, lambda_mid, valNELE, matPLEX, matEATT, matELOC);

%% Vorticity

% Typically found at the two vertices that split two HDVEs
% vnuma is local vertices 1 and 2 (columns) for HDVE 1
% vnumb is local vertices 1 and 2 for HDVE 2
[ta,tb,~] = find(matDVE(matEATT(idx,1),:,1) == repmat(matELST(idx,1),1,3));
[~,rb,] = sort(ta);
vnuma(:,1) = tb(rb); % Sorting it according to row, cause find returns them jumbled up

[ta,tb,~] = find(matDVE(matEATT(idx,1),:,1) == repmat(matELST(idx,2),1,3));
[~,rb,] = sort(ta);
vnuma(:,2) = tb(rb);

[ta,tb,~] = find(matDVE(matEATT(idx,2),:,1) == repmat(matELST(idx,1),1,3));
[~,rb,] = sort(ta);
vnumb(:,1) = tb(rb);

[ta,tb,~] = find(matDVE(matEATT(idx,2),:,1) == repmat(matELST(idx,2),1,3));
[~,rb,] = sort(ta);
vnumb(:,2) = tb(rb);

% First vertex
[vort_e1, vort_x1] = fcnDVORT(idx, vnuma(:,1), vnumb(:,1), nedg, lambda_vert, valNELE, matPLEX, matEATT, matELOC, matALIGN);

% Second vertex
[vort_e2, vort_x2] = fcnDVORT(idx, vnuma(:,2), vnumb(:,2), nedg, lambda_vert, valNELE, matPLEX, matEATT, matELOC, matALIGN);

%% Irrotationality

dvenum = [1:valNELE]';

% A1 + B1 = 0
ddgamma = repmat([1 0 -1 0 0],valNELE,1);

rows = reshape([repmat([1:valNELE]',1,5)]',[],1);
col3 = reshape([repmat((dvenum.*5)-4,1,5) + repmat([0:4],valNELE,1)]',[],1);

irrot = zeros(valNELE, valNELE*5);
irrot(sub2ind(size(irrot),rows,col3)) = reshape(ddgamma',[],1);

%% Circulation equations at leading edge and wing tip
% For lifting surface analysis
% Circulation is set to zero at the leading edge, and wing tips
% These are found by looking at the free edges that are NOT symmetry or trailing edge
% Evaluated at the mid-point of each edge which is used by only 1 HDVE (and not at the trailing edge)
circ_tip = [];
if strcmp(strATYPE,'LS') == 1
    idx = ~all(matEATT,2); % All edges that are attached to only 1 HDVE
    idx(vecTE) = 0;
    idx(vecLE) = 0;
    
    nedg = length(matEATT(idx,1));
    
    circ_tip = fcnDCIRCTIP(idx, nedg, lambda_mid, valNELE, matPLEX, matEATT, matELOC);   
end

%% Kinematic conditions at vertices
% Flow tangency is to be enforced at all control points on the surface HDVEs

% In the D-Matrix, dot (a1,a2,b1,b2,c3) of our influencing HDVE with the normal of the point we are influencing on

% Points we are influencing
% fpg = [VLST; CENTER];
fpg = [matCENTER];

% List of DVEs we are influencing from (one for each of the above fieldpoints)
len = length(fpg(:,1));
dvenum = reshape(repmat(1:valNELE,len,1),[],1);
dvetype = ones(size(dvenum));

fpg = repmat(fpg,valNELE,1);

[a1, a2, b1, b2, c3] = fcnHDVEIND(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX, dvetype);

% List of normals we are to dot the above with
% normals = [VNORM; DVECT(:,:,3)];
normals = [matDVECT(:,:,3)];
normals = repmat(normals,valNELE,1); % Repeated so we can dot all at once

% Dotting a1, a2, b1, b2, c3 with the normals of the field points
temp60 = [dot(a1,normals,2) dot(a2,normals,2) dot(b1,normals,2) dot(b2,normals,2) dot(c3, normals,2)];

% Reshaping and inserting into the bottom of the D-Matrix
rows = [1:len]';

king_kong = zeros(len, valNELE*5);
king_kong(rows,:) = reshape(permute(reshape(temp60',5,[],valNELE),[2 1 3]),[],5*valNELE,1);

%% Piecing together D-matrix

D = [circ; vort_e1; vort_x1; vort_e2; vort_x2; irrot; circ_tip; king_kong];

end

