function [D] = fcnDWING9(strATYPE, matEATT, matPLEX, valNELE, matELOC, matELST, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecLEDVE, vecTEDVE, vecSYM, matROTANG, vecSPANDIR, matKINCON_P, matKINCON_DVE)

%% Circulation equations between elements
% Evaluated at the mid-point of each edge which splits two HDVEs
idx = all(matEATT,2); % All edges that split 2 DVEs

if ~isempty(vecTE)
   idx(vecTE) = false; 
end

vnum_a = matVLST(matELST(idx,1),:);
vnum_b = matVLST(matELST(idx,2),:);
vnum_mid = (vnum_a + vnum_b)./2;

%% Circulation at edge corner and midpoints
dvenum = [matEATT(idx,1) matEATT(idx,2)];
circ = [fcnDCIRC(repmat(vnum_a,1,1,2), dvenum, valNELE, matROTANG, matCENTER); ...
        fcnDCIRC(repmat(vnum_mid,1,1,2), dvenum, valNELE, matROTANG, matCENTER); ...
        fcnDCIRC(repmat(vnum_b,1,1,2), dvenum, valNELE, matROTANG, matCENTER)];

%% Vorticity along edge between elements
% Unit vector in local ref frame (a for HDVE1, b for HDVE2) from local vertex to local vertex on the edge that forms the border between the two
vort = [fcnDVORTEDGE(repmat(vnum_a,1,1,2), dvenum, valNELE, matROTANG, matCENTER); ...
        fcnDVORTEDGE(repmat(vnum_mid,1,1,2), dvenum, valNELE, matROTANG, matCENTER); ...
        fcnDVORTEDGE(repmat(vnum_b,1,1,2), dvenum, valNELE, matROTANG, matCENTER)];

%% Circulation equations at wing tip (and LE?)
% For lifting surface analysis
% Circulation is set to zero at the wing tips
% These are found by looking at the free edges that are NOT symmetry or trailing edge
% Evaluated at the mid-point of each edge which is used by only 1 HDVE (and not at the trailing edge)
circ_tip = [];
vort_tip = [];
vort_te = [];

if strcmpi(strATYPE{2},'THIN') == 1
    pts(:,:,1) = matVLST(matELST(vecLE,1),:);
    pts(:,:,2) = matVLST(matELST(vecLE,2),:);
    pts(:,:,3) = (pts(:,:,1) + pts(:,:,2))./2;
    
    circ_tip = [fcnDCIRC2(pts(:,:,1), vecLEDVE, valNELE, matROTANG, matCENTER); ...
                fcnDCIRC2(pts(:,:,2), vecLEDVE, valNELE, matROTANG, matCENTER); ...
                fcnDCIRC2(pts(:,:,3), vecLEDVE, valNELE, matROTANG, matCENTER)];    
 
    pts = [];
    pts(:,:,1) = matVLST(matELST(vecTE,1),:);
    pts(:,:,2) = matVLST(matELST(vecTE,2),:);
    pts(:,:,3) = (pts(:,:,1) + pts(:,:,2))./2;
    
    te_dir = pts(:,:,2) - pts(:,:,1);
    te_dir = te_dir./(sqrt(sum(te_dir.^2,2)));
    
%     vort_te = [te_dir(:,2).*fcnDVORT2(pts(:,:,1), vecTEDVE, valNELE, matCENTER, matROTANG, 'A') + te_dir(:,1).*fcnDVORT2(pts(:,:,1), vecTEDVE, valNELE, matCENTER, matROTANG, 'B'); ...
%     te_dir(:,2).*fcnDVORT2(pts(:,:,2), vecTEDVE, valNELE, matCENTER, matROTANG, 'A') + te_dir(:,1).*fcnDVORT2(pts(:,:,2), vecTEDVE, valNELE, matCENTER, matROTANG, 'B'); ...
%     te_dir(:,2).*fcnDVORT2(pts(:,:,3), vecTEDVE, valNELE, matCENTER, matROTANG, 'A') + te_dir(:,1).*fcnDVORT2(pts(:,:,3), vecTEDVE, valNELE, matCENTER, matROTANG, 'B')];
    
    vort_te = [fcnDVORT2(pts(:,:,1), vecTEDVE, valNELE, matCENTER, matROTANG, 'A'); ...
        fcnDVORT2(pts(:,:,2), vecTEDVE, valNELE, matCENTER, matROTANG, 'A'); ...
        fcnDVORT2(pts(:,:,3), vecTEDVE, valNELE, matCENTER, matROTANG, 'A')];
end
            
if strcmpi(strATYPE{1},'3D') == 1    
    idx = ~all(matEATT,2);
    idx(vecTE) = 0;
    idx(vecLE) = 0;
    pts2(:,:,1) = matVLST(matELST(idx,1),:);
    pts2(:,:,2) = matVLST(matELST(idx,2),:);
    pts2(:,:,3) = (pts2(:,:,1) + pts2(:,:,2))./2;
    tipdves = nonzeros(matEATT(idx,:));
    
    circ_tip = [circ_tip; 
                fcnDCIRC2(pts2(:,:,1), tipdves, valNELE, matROTANG, matCENTER); ...
                fcnDCIRC2(pts2(:,:,2), tipdves, valNELE, matROTANG, matCENTER); ...
                fcnDCIRC2(pts2(:,:,3), tipdves, valNELE, matROTANG, matCENTER)];
elseif strcmpi(strATYPE{1},'2D') == 1
    dvenum = [1:valNELE]';    
    vort_tip = fcnDVORT1(dvenum, valNELE,'B');
end

%% Trailing edge vorticity

%% Kinematic conditions at vertices
% Flow tangency is to be enforced at all control points on the surface HDVEs
% In the D-Matrix, dot (a1,a2,b1,b2,c3) of our influencing HDVE with the normal of the point we are influencing on

% List of DVEs we are influencing from (one for each of the above fieldpoints)
len = length(matKINCON_P(:,1));
dvenum = reshape(repmat(1:valNELE,len,1),[],1);
dvetype = ones(size(dvenum));

fpg = repmat(matKINCON_P,valNELE,1);

[infl_glob] = fcnHDVEINDGLOB(dvenum, dvetype, fpg, matPLEX, matROTANG, matCENTER);

normals = repmat(matDVECT(matKINCON_DVE,:,3),valNELE,1); % Repeated so we can dot all at once

% Dotting a1, a2, b1, b2, c3 with the normals of the field points
temp60 = [dot(permute(infl_glob(:,1,:),[3 1 2]),normals,2) dot(permute(infl_glob(:,2,:),[3 1 2]),normals,2) dot(permute(infl_glob(:,3,:),[3 1 2]),normals,2) dot(permute(infl_glob(:,4,:),[3 1 2]),normals,2) dot(permute(infl_glob(:,5,:),[3 1 2]),normals,2)];

% Reshaping and inserting into the bottom of the D-Matrix
rows = [1:len]';

king_kong = zeros(len, valNELE*5);
king_kong(rows,:) = reshape(permute(reshape(temp60',5,[],valNELE),[2 1 3]),[],5*valNELE,1);

%% Piecing together D-matrix
D = [circ; vort; vort_tip; vort_te; circ_tip; king_kong];
% D = [circ; vort; vort_te; king_kong];
end

