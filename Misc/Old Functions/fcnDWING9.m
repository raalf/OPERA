function [D] = fcnDWING9(strATYPE, matEATT, matPLEX, valNELE, matELOC, matELST, matALIGN, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecSYM, matVATT, vecTEDVE, vecSPANDIR, matROTANG, matVNORM)

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

% Making D-matrix a little long, in case we are over-constrained
D = zeros(valNELE*9, valNELE*6);

%% Circulation equations between elements

% Evaluated at the mid-point of each edge which splits two HDVEs
idx = all(matEATT,2); % All edges that split 2 DVEs
nedg = length(matEATT(idx,1));

dve1 = matEATT(idx,1);
dve2 = matEATT(idx,2);

step = (matVNORM(matELST(idx,1),:) + matVNORM(matELST(idx,2),:))./2;
fpg = (matVLST(matELST(idx,1),:) + matVLST(matELST(idx,2),:))./2; 
fpg = fpg + 0.1.*step;

dinf1 = fcnDINF(idx, nedg, dve1, dve2, fpg, matDVE, matDVECT, matVLST, matPLEX, ones(length(dve1),1), matROTANG, matEATT, valNELE);


%% Vorticity along edge between elements
step = matVNORM(matELST(idx,1),:);
fpg = matVLST(matELST(idx,1),:);
fpg = fpg + -0.1.*step;
dinf2 = fcnDINF(idx, nedg, dve1, dve2, fpg, matDVE, matDVECT, matVLST, matPLEX, ones(length(dve1),1), matROTANG, matEATT, valNELE);

step = matVNORM(matELST(idx,2),:);
fpg = matVLST(matELST(idx,2),:);
fpg = fpg + -0.1.*step;
dinf3 = fcnDINF(idx, nedg, dve1, dve2, fpg, matDVE, matDVECT, matVLST, matPLEX, ones(length(dve1),1), matROTANG, matEATT, valNELE);

%% Circulation equations at wing tip (and LE?)
% For lifting surface analysis
% Circulation is set to zero at the wing tips
% These are found by looking at the free edges that are NOT symmetry or trailing edge
% Evaluated at the mid-point of each edge which is used by only 1 HDVE (and not at the trailing edge)
circ_tip = [];
if strcmp(strATYPE,'THIN') == 1
    idx = ~all(matEATT,2); % All edges that are attached to only 1 HDVE
    idx(vecTE) = 0;
    idx(vecLE) = 0;
    
    nedg = length(matEATT(idx,1));
    
    % Finding the two local vertex numbers of the edge endpoints
    vnumc = nonzeros(matELOC(idx,:));
    vnumc = repmat(vnumc,1,2);
    vnumc(:,2) = vnumc(:,2) + 1;
    vnumc(vnumc(:,2) == 4,2) = 1;
    
    circ_tip = fcnDCIRCTIP(idx, nedg, lambda_vert, lambda_mid, valNELE, matPLEX, matEATT, matELOC, vnumc);  
end

%% Kinematic conditions at vertices
% Flow tangency is to be enforced at all control points on the surface HDVEs

% In the D-Matrix, dot (a1,a2,b1,b2,c3) of our influencing HDVE with the normal of the point we are influencing on

% Points we are influencing
fpg = matCENTER;
normals = matDVECT(:,:,3);

% List of DVEs we are influencing from (one for each of the above fieldpoints)
len = length(fpg(:,1));
dvenum = reshape(repmat(1:valNELE,len,1),[],1);
dvetype = ones(size(dvenum));

fpg = repmat(fpg,valNELE,1);

[a1, a2, a3, b1, b2, b3] = fcnHDVEIND(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX, dvetype, matROTANG);

normals = repmat(normals,valNELE,1); % Repeated so we can dot all at once

% Dotting a1, a2, b1, b2, c3 with the normals of the field points
temp60 = [dot(a1,normals,2) dot(a2,normals,2) dot(a3,normals,2) dot(b1,normals,2) dot(b2,normals,2) dot(b3,normals,2)];

% Reshaping and inserting into the bottom of the D-Matrix
rows = [1:len]';

king_kong = zeros(len, valNELE*6);
king_kong(rows,:) = reshape(permute(reshape(temp60',6,[],valNELE),[2 1 3]),[],6*valNELE,1);

%% Piecing together D-matrix
D = [dinf1; dinf2; dinf3; circ_tip; king_kong];

end

