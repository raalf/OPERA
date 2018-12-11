function [D] = fcnDWING9(strATYPE, matEATT, matPLEX, valNELE, matELOC, matELST, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecSYM, matROTANG, vecSPANDIR, vecTEDVE, matCONTROL)

%% Circulation equations between elements
% Evaluated at the mid-point of each edge which splits two HDVEs
idx = all(matEATT,2); % All edges that split 2 DVEs
vnum_a = matVLST(matELST(idx,1),:);
vnum_b = matVLST(matELST(idx,2),:);
vnum_mid = (vnum_a + vnum_b)./2;

%% Circulation at edge corner and midpoints
dvenum = [matEATT(idx,1) matEATT(idx,2)];
circ = [fcnDCIRC(idx, repmat(vnum_a,1,1,2), dvenum, valNELE, matROTANG, matEATT, matCONTROL); ...
        fcnDCIRC(idx, repmat(vnum_mid,1,1,2), dvenum, valNELE, matROTANG, matEATT, matCONTROL); ...
        fcnDCIRC(idx, repmat(vnum_b,1,1,2), dvenum, valNELE, matROTANG, matEATT, matCONTROL)];

%% Vorticity along edge between elements
% Unit vector in local ref frame (a for HDVE1, b for HDVE2) from local vertex to local vertex on the edge that forms the border between the two

xi_in_ngbr = fcnGLOBSTAR(matDVECT(dvenum(:,1),:,1), matROTANG(dvenum(:,2),:));
xi_in_ngbr = xi_in_ngbr - (dot(xi_in_ngbr, matDVECT(dvenum(:,2),:,3), 2)).*xi_in_ngbr;
xi_in_ngbr = xi_in_ngbr./(sqrt(sum(xi_in_ngbr.^2, 2)));
theta = acos(dot(xi_in_ngbr, repmat([1 0 0], size(xi_in_ngbr,1), 1), 2));

vort = [];
% vort = [fcnDVORTEDGE(idx, repmat(vnum_a,1,1,2), theta, dvenum, valNELE, matROTANG, matEATT, matCENTER); ...
%         fcnDVORTEDGE(idx, repmat(vnum_mid,1,1,2), theta, dvenum, valNELE, matROTANG, matEATT, matCENTER); ...
%         fcnDVORTEDGE(idx, repmat(vnum_b,1,1,2), theta, dvenum, valNELE, matROTANG, matEATT, matCENTER)];

%% Circulation equations at wing tip (and LE?)
% For lifting surface analysis
% Circulation is set to zero at the wing tips
% These are found by looking at the free edges that are NOT symmetry or trailing edge
% Evaluated at the mid-point of each edge which is used by only 1 HDVE (and not at the trailing edge)
circ_tip = [];
vort_tip = [];

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

    idx = ~all(matEATT,2); % All edges that are attached to only 1 HDVE
    idx(vecTE) = 0;
    idx(vecLE) = 0;
    
    dvenum = [1:valNELE]';
    theta = acos(dot(matDVECT(dvenum,:,2), vecSPANDIR(dvenum,:),2));
        
    pts(:,:,1) = matCENTER(dvenum,:) + vecSPANDIR(dvenum,:);
    pts(:,:,2) = matCENTER(dvenum,:) + matDVECT(dvenum,:,2);
    pts(:,:,3) = matCENTER(dvenum,:) + matDVECT(dvenum,:,1);
%     vort_tip = [fcnDVORTTE(idx, pts(:,:,1), theta, dvenum, valNELE, matROTANG, matEATT, matCENTER); ...
%                 fcnDVORTTE(idx, pts(:,:,2), theta, dvenum, valNELE, matROTANG, matEATT, matCENTER); ...
%                 fcnDVORTTE(idx, pts(:,:,3), theta, dvenum, valNELE, matROTANG, matEATT, matCENTER)];
    
    vort_tip = [fcnDVORTTE(idx, pts(:,:,2), theta, dvenum, valNELE, matROTANG, matEATT, matCENTER); ...
                fcnDVORTTE(idx, pts(:,:,3), theta, dvenum, valNELE, matROTANG, matEATT, matCENTER)];

%     dvenum = repmat(nonzeros(matEATT(idx,:)), 1, 2);
%     vnum_a = (matVLST(matELST(idx,1),:) + matVLST(matELST(idx,2),:))./2;  
% %     vnum_b = vnum_a + (dot(vecSPANDIR(dvenum(:,1),:), matCENTER(dvenum(:,1),:) - vnum_a,2)./(sqrt(sum(vecSPANDIR(dvenum(:,1),:).^2,2))).^2).*vecSPANDIR(dvenum(:,1),:);
%     vnum_b = vnum_a + vecSPANDIR(dvenum(:,1),:);
%     vnum_mid = (vnum_a + vnum_b)./2;
%        
%     vort_tip = [fcnDCIRC(idx, cat(3, vnum_a, vnum_mid), dvenum, valNELE, matROTANG, matEATT, matCONTROL); ...
%                 fcnDCIRC(idx, cat(3, vnum_mid, vnum_b), dvenum, valNELE, matROTANG, matEATT, matCONTROL)];     
end

%% Trailing edge vorticity

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
temp60 = [dot(permute(infl_glob(:,1,:),[3 1 2]),normals,2) dot(permute(infl_glob(:,2,:),[3 1 2]),normals,2) dot(permute(infl_glob(:,3,:),[3 1 2]),normals,2) dot(permute(infl_glob(:,4,:),[3 1 2]),normals,2) dot(permute(infl_glob(:,5,:),[3 1 2]),normals,2)];

% Reshaping and inserting into the bottom of the D-Matrix
rows = [1:len]';

king_kong = zeros(len, valNELE*5);
king_kong(rows,:) = reshape(permute(reshape(temp60',5,[],valNELE),[2 1 3]),[],5*valNELE,1);

%% Piecing together D-matrix
D = [circ; vort; vort_tip; circ_tip; king_kong];

end

