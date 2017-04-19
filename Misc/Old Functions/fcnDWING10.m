function [D] = fcnDWING10(strATYPE, matEATT, matPLEX, valNELE, matELOC, matELST, matALIGN, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecSYM, matVATT, vecTEDVE, vecSPANDIR, matROTANG)

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


%% Circulation equations between elements

% Evaluated at the mid-point of each edge which splits two HDVEs
idx = all(matEATT,2); % All edges that split 2 DVEs
nedg = length(matEATT(idx,1));

circ = fcnDCIRC(idx, nedg, lambda_mid, valNELE, matPLEX, matEATT, matELOC);

%% Vorticity

% Setting vorticity in the eta and xsi directions equal between neighbouring HDVEs, at edge mid points
[vort_e, vort_x] = fcnDVORT2(idx, nedg, lambda_mid, valNELE, matPLEX, matEATT, matELOC, matALIGN);

%% Setting curvature equal at edge mid points

[curv] = fcnDCURV(idx, nedg, lambda_mid, valNELE, matPLEX, matEATT, matELOC);

%% Irrotationality

dvenum = [1:valNELE]';

% A1 + B1 = 0
ddgamma = repmat([1 0 -1 0 0],valNELE,1);

rows = reshape([repmat([1:valNELE]',1,5)]',[],1);
col3 = reshape([repmat((dvenum.*5)-4,1,5) + repmat([0:4],valNELE,1)]',[],1);

irrot = zeros(valNELE, valNELE*5);
irrot(sub2ind(size(irrot),rows,col3)) = reshape(ddgamma',[],1);

%% Circulation equations at wing tip
% For lifting surface analysis
% Circulation is set to zero at the wing tips
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

%% Setting TE Spanwise Vorticity 

vort_te = [];

if ~isempty(vecTEDVE)
    len = length(vecTEDVE);

    %(x,y) of all three vertices of HDVEs at trailing edge of the wing and the corresponding leading edge of the post te wake HDVE
    x1 = reshape(matPLEX(1,1,vecTEDVE),len,1);
    x2 = reshape(matPLEX(2,1,vecTEDVE),len,1);
    x3 = reshape(matPLEX(3,1,vecTEDVE),len,1);
    y1 = reshape(matPLEX(1,2,vecTEDVE),len,1);
    y2 = reshape(matPLEX(2,2,vecTEDVE),len,1);
    y3 = reshape(matPLEX(3,2,vecTEDVE),len,1);

    lmb1 = reshape(lambda_mid([nonzeros(matELOC(vecTE,:))],1),len,1);
    lmb2 = reshape(lambda_mid([nonzeros(matELOC(vecTE,:))],2),len,1);
    lmb3 = reshape(lambda_mid([nonzeros(matELOC(vecTE,:))],3),len,1);

    a2 = ones(len,1);
    a1 = (lmb1.*x1+lmb2.*x2+lmb3.*x3);
    b2 = ones(len,1);
    b1 = (lmb1.*y1+lmb2.*y2+lmb3.*y3);

    c3 = zeros(len,1);

    dgammate = [a1(:,1).*vecSPANDIR(:,2), a2(:,1).*vecSPANDIR(:,2), b1(:,1).*vecSPANDIR(:,1), b2(:,1).*vecSPANDIR(:,1), c3(:,1)];

    % Row indices of the rows where vorticity equations will go
    rows = reshape([repmat([1:len]',1,5)]',[],1);

    % Column indices for each circulation equation, col# = (DVE*6)-5 as each DVE gets a 6 column group
    col1 = reshape([repmat([(vecTEDVE.*5)-4],1,5)+repmat([0:4],len,1)]',[],1);

    vort_te = zeros(len, valNELE*5);
    vort_te(sub2ind(size(vort_te),rows,col1)) = reshape(dgammate,[],1);
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

[a1, a2, b1, b2, c3] = fcnHDVEIND(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX, dvetype, matROTANG);

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

D = [circ; vort_e; vort_x; curv; irrot; circ_tip; king_kong];
end

