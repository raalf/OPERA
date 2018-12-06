function [D] = fcnDWING9(strATYPE, matEATT, matPLEX, valNELE, matELOC, matELST, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecSYM, matROTANG, vecSPANDIR, vecTEDVE, matCONTROL)

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

%% Circulation
circ = [];
circ1 = [];
circ2 = [];

circ = fcnDCIRC(idx, nedg, lambda_mid, valNELE, matPLEX, matEATT, matELOC);
circ1 = fcnDCIRC2(idx, vnuma(:,1), vnumb(:,1), nedg, lambda_vert, valNELE, matPLEX, matEATT, matELOC);
circ2 = fcnDCIRC2(idx, vnuma(:,2), vnumb(:,2), nedg, lambda_vert, valNELE, matPLEX, matEATT, matELOC);

%% Vorticity along edge between elements
% Unit vector in local ref frame (a for HDVE1, b for HDVE2) from local vertex to local vertex on the edge that forms the border between the two
vort1 = [];
vort2 = [];
vort3 = [];
vort4 = [];

temp = reshape(permute(matPLEX,[1 3 2]),[],3,1);

e1vec = temp(matEATT(idx,1).*3 + vnuma(:,2) - 3,:) - temp(matEATT(idx,1).*3 + vnuma(:,1) - 3,:);
e2vec = temp(matEATT(idx,2).*3 + vnumb(:,2) - 3,:) - temp(matEATT(idx,2).*3 + vnumb(:,1) - 3,:);
e1vec = e1vec(:,1:2)./repmat(sqrt(e1vec(:,1).^2 + e1vec(:,2).^2),1,2);
e2vec = e2vec(:,1:2)./repmat(sqrt(e2vec(:,1).^2 + e2vec(:,2).^2),1,2);

vort1 = fcnDVORTEDGE(idx, vnuma(:,1), vnumb(:,1), nedg, lambda_vert, valNELE, matPLEX, matEATT, e1vec, e2vec);
vort2 = fcnDVORTEDGE(idx, vnuma(:,2), vnumb(:,2), nedg, lambda_vert, valNELE, matPLEX, matEATT, e1vec, e2vec);

vort3 = fcnDVORTEDGE(idx, vnuma(:,1), vnumb(:,1), nedg, lambda_vert, valNELE, matPLEX, matEATT, [-e1vec(:,2) e1vec(:,1)], [-e2vec(:,2) e2vec(:,1)]);
vort4 = fcnDVORTEDGE(idx, vnuma(:,2), vnumb(:,2), nedg, lambda_vert, valNELE, matPLEX, matEATT, [-e1vec(:,2) e1vec(:,1)], [-e2vec(:,2) e2vec(:,1)]);

%% Circulation equations at wing tip (and LE?)
% For lifting surface analysis
% Circulation is set to zero at the wing tips
% These are found by looking at the free edges that are NOT symmetry or trailing edge
% Evaluated at the mid-point of each edge which is used by only 1 HDVE (and not at the trailing edge)
circ_tip = [];
vort_tip1 = [];
vort_tip2 = [];

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

elseif strcmp(strATYPE,'2D') == 1

    vnuma = [];
    
    idx = ~all(matEATT,2); % All edges that are attached to only 1 HDVE
    idx(vecTE) = 0;
    idx(vecLE) = 0;
    
    nedg = length(matEATT(idx,1));
    
    [ta,tb,~] = find(matDVE(nonzeros(matEATT(idx,:)),:) == repmat(matELST(idx,1),1,3));
    [~,rb,] = sort(ta);
    vnuma(:,1) = tb(rb); 

    [ta,tb,~] = find(matDVE(nonzeros(matEATT(idx,:)),:) == repmat(matELST(idx,2),1,3));
    [~,rb,] = sort(ta);
    vnuma(:,2) = tb(rb);

    e1vec = temp(nonzeros(matEATT(idx,:)).*3 + vnuma(:,2) - 3,:) - temp(nonzeros(matEATT(idx,:)).*3 + vnuma(:,1) - 3,:);
    e1vec = e1vec./sqrt(sum(e1vec.^2,2));
    e1vec = [-e1vec(:,2), e1vec(:,1), e1vec(:,3)];
    vort_tip1 = fcnDVORTTE(idx, vnuma(:,1), sum(idx), lambda_vert, valNELE, matPLEX, matEATT, e1vec);
    vort_tip2 = fcnDVORTTE(idx, vnuma(:,2), sum(idx), lambda_vert, valNELE, matPLEX, matEATT, e1vec);
    
end

%% Trailing edge vorticity
vort5 = [];
vort6 = [];
vnuma = [];
vnumb = [];

% % [ta,tb,~] = find(matDVE(vecTEDVE,:) == repmat(matELST(vecTE,1),1,3));
% % [~,rb,] = sort(ta);
% % vnuma(:,1) = tb(rb); 
% % 
% % [ta,tb,~] = find(matDVE(vecTEDVE,:) == repmat(matELST(vecTE,2),1,3));
% % [~,rb,] = sort(ta);
% % vnuma(:,2) = tb(rb);
% % 
% % % e1vec = temp(nonzeros(matEATT(vecTE,:)).*3 + vnuma(:,2) - 3,:) - temp(nonzeros(matEATT(vecTE,:)).*3 + vnuma(:,1) - 3,:);
% % e1vec = vecSPANDIR;
% % % vort5 = fcnDVORTTE(vecTE, vnuma(:,1), length(vecTE), lambda_vert, valNELE, matPLEX, matEATT, e1vec);
% % % vort6 = fcnDVORTTE(vecTE, vnuma(:,2), length(vecTE), lambda_vert, valNELE, matPLEX, matEATT, e1vec);

%% Kinematic conditions at vertices
% Flow tangency is to be enforced at all control points on the surface HDVEs
% In the D-Matrix, dot (a1,a2,b1,b2,c3) of our influencing HDVE with the normal of the point we are influencing on

% Points we are influencing
fpg = matCENTER;
normals = matDVECT(:,:,3);

% fpg = matCENTER + (matDVECT(:,:,3)./10000).*-1;

% List of DVEs we are influencing from (one for each of the above fieldpoints)
len = length(fpg(:,1));
dvenum = reshape(repmat(1:valNELE,len,1),[],1);
dvetype = ones(size(dvenum));

fpg = repmat(fpg,valNELE,1);

[infl_glob] = fcnHDVEIND(dvenum, dvetype, fpg, matPLEX, matROTANG, matCONTROL);

normals = repmat(normals,valNELE,1); % Repeated so we can dot all at once

% Dotting a1, a2, b1, b2, c3 with the normals of the field points
temp60 = [dot(permute(infl_glob(:,1,:),[3 1 2]),normals,2) dot(permute(infl_glob(:,2,:),[3 1 2]),normals,2) dot(permute(infl_glob(:,3,:),[3 1 2]),normals,2) dot(permute(infl_glob(:,4,:),[3 1 2]),normals,2) dot(permute(infl_glob(:,5,:),[3 1 2]),normals,2) dot(permute(infl_glob(:,6,:),[3 1 2]),normals,2)];

% Reshaping and inserting into the bottom of the D-Matrix
rows = [1:len]';

king_kong = zeros(len, valNELE*6);
king_kong(rows,:) = reshape(permute(reshape(temp60',6,[],valNELE),[2 1 3]),[],6*valNELE,1);

%% Piecing together D-matrix
D = [circ; circ1; circ2; vort1; vort2; vort3; vort4; vort5; vort6; vort_tip1; vort_tip2; circ_tip; king_kong];
% D = [circ; circ1; circ2; vort1; vort2; vort3; vort4; vort5; vort6; vort_tip1; vort_tip2; circ_tip; king_kong];

end

