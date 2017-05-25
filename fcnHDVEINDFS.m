function [a1, a2, a3, b1, b2, b3] = fcnHDVEINDFS(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX, dvetype, matROTANG, matVSCOMB)
% dvetype - 1 for surface, 2 for wake, 3 for semi-infinite wake

dvenum = reshape(dvenum, [], 1, 1); % Ensuring dvenum is a column vector

len = length(dvenum);
dbl_eps = 1e-14;

fp = fcnGLOBSTAR(fpg - matVLST(matDVE(dvenum,1,1),:), matROTANG(dvenum,1), matROTANG(dvenum,2), matROTANG(dvenum,3));

endpoints = zeros(len*6, 3, 2); % Five function calls per DVE
phi = zeros(len*6,1);
yaw = zeros(len*6,1);

k = zeros(len*6,1) + 0.1; % Until this is figured out for HDVE, k = 1

idx = [1:6:len*6]';

%% Finding velocities induced by each vortex sheet
% Always input left point, right point of vortex sheet LE

% Right to left. Sheet LE is on right, sheet extends to infinity on the left.
% First edge
endpoints(idx,:,1:2) = reshape(permute(matPLEX(1:2,:,dvenum),[3 2 1]), [],3,2);
yaw(idx) = pi;
phi(idx) = atan((endpoints(idx,1,2)-endpoints(idx,1,1))./(endpoints(idx,2,2)-endpoints(idx,2,1)));

% Second edge
endpoints(idx+1,:,1:2) = reshape(permute(matPLEX(2:3,:,dvenum),[3 2 1]), [],3,2);
yaw(idx+1) = pi;
phi(idx+1) = atan((endpoints(idx+1,1,2)-endpoints(idx+1,1,1))./(endpoints(idx+1,2,2)-endpoints(idx+1,2,1)));

% Third edge
endpoints(idx+2,:,1:2) = reshape(permute(matPLEX(1:2:3,:,dvenum),[3 2 1]), [],3,2);
yaw(idx+2) = pi;
phi(idx+2) = atan((endpoints(idx+2,1,2)-endpoints(idx+2,1,1))./(endpoints(idx+2,2,2)-endpoints(idx+2,2,1)));

% Up to down. Sheet LE is on the bottom, sheet extends to infinity upwards.
% First edge
endpoints(idx+3,:,1:2) = reshape(permute(matPLEX(1:2,:,dvenum),[3 2 1]), [],3,2);
yaw(idx+2) = pi/2;
phi(idx+2) = -atan((endpoints(idx+2,2,2)-endpoints(idx+2,2,1))./(endpoints(idx+2,1,2)-endpoints(idx+2,1,1)));

% Second edge
endpoints(idx+4,:,1:2) = reshape(permute(matPLEX(2:3,:,dvenum),[3 2 1]), [],3,2);
yaw(idx+3) = pi/2;
phi(idx+3) = -atan((endpoints(idx+3,2,2)-endpoints(idx+3,2,1))./(endpoints(idx+3,1,2)-endpoints(idx+3,1,1)));

% Third edge
endpoints(idx+5,:,1:2) = reshape(permute(matPLEX(1:2:3,:,dvenum),[3 2 1]), [],3,2);
yaw(idx+4) = pi/2;
phi(idx+4) = -atan((endpoints(idx+5,2,2)-endpoints(idx+5,2,1))./(endpoints(idx+5,1,2)-endpoints(idx+5,1,1))); % Leading edge of Edge 3 is always on eta-axis

%% Running VSIND

fpl = reshape(repmat(fp,1,6,1)',3,[],1)';

[al, bl, cl] = fcnVSIND(endpoints, phi, yaw, fpl, k);

a1l = zeros(len,3);
a2l = zeros(len,3);
a3l = zeros(len,3);

%% Summing the velocities of all five sheets

% RIGHT-TO-LEFT
b1l = repmat(matVSCOMB(dvenum,1,1),1,3).*cl(idx,:) + repmat(matVSCOMB(dvenum,2,1),1,3).*cl(idx+1,:) + repmat(matVSCOMB(dvenum,3,1),1,3).*cl(idx+2,:);
b2l = repmat(matVSCOMB(dvenum,1,1),1,3).*bl(idx,:) + repmat(matVSCOMB(dvenum,2,1),1,3).*bl(idx+1,:) + repmat(matVSCOMB(dvenum,3,1),1,3).*bl(idx+2,:);
b3l = repmat(matVSCOMB(dvenum,1,1),1,3).*al(idx,:) + repmat(matVSCOMB(dvenum,2,1),1,3).*al(idx+1,:) + repmat(matVSCOMB(dvenum,3,1),1,3).*al(idx+2,:);

% DOWN-TO-UP (negative because matVSCOMB assumes sheet goes left to right, to positive eta infinity)
a1l = -repmat(matVSCOMB(dvenum,1,2),1,3).*cl(idx+3,:) + -repmat(matVSCOMB(dvenum,2,2),1,3).*cl(idx+4,:) + -repmat(matVSCOMB(dvenum,3,2),1,3).*cl(idx+5,:);
a2l = -repmat(matVSCOMB(dvenum,1,2),1,3).*bl(idx+3,:) + -repmat(matVSCOMB(dvenum,2,2),1,3).*bl(idx+4,:) + -repmat(matVSCOMB(dvenum,3,2),1,3).*bl(idx+5,:);
a3l = -repmat(matVSCOMB(dvenum,1,2),1,3).*al(idx+3,:) + -repmat(matVSCOMB(dvenum,2,2),1,3).*al(idx+4,:) + -repmat(matVSCOMB(dvenum,3,2),1,3).*al(idx+5,:);


% % Subtracting first vortex sheet from second in the right-to-left direction
% b1l = cl(idx+1,:) - cl(idx,:);
% b2l = bl(idx+1,:) - bl(idx,:);
% b3l = al(idx+1,:) - al(idx,:);
% 
% % % Leaving the semi-infites on the oldest wake row
% b1l(dvetype == 3,:) = cl(idx(dvetype == 3) + 1,:);
% b2l(dvetype == 3,:) = bl(idx(dvetype == 3) + 1,:);
% b3l(dvetype == 3,:) = al(idx(dvetype == 3) + 1,:);
% 
% % Subtracting the three vortex sheets in the up-to-down direction based on the shape of the triangle
% 
% idx_a = endpoints(idx,1,2) < 0; % HDVEs with obtuse angle at Vertex 1 (eta < 0)
% idx_a1 = idx(idx_a.*idx ~= 0); % Index of first edge of HDVEs for obtuse angle at Vertex 1
% a1l(idx_a,:) = cl(idx_a1+2,:) - cl(idx_a1+3,:) + cl(idx_a1+4,:);
% a2l(idx_a,:) = bl(idx_a1+2,:) - bl(idx_a1+3,:) + bl(idx_a1+4,:);
% a3l(idx_a,:) = al(idx_a1+2,:) - al(idx_a1+3,:) + al(idx_a1+4,:);
% 
% idx_b = endpoints(idx,1,2) >= 0 & endpoints(idx+1,1,1) <= endpoints(idx+1,1,2); % HDVEs with acute (or right) angle at Vertex 1 & Vertex 2
% idx_b1 = idx(idx_b.*idx ~= 0);
% a1l(idx_b,:) = -cl(idx_b1+2,:) - cl(idx_b1+3,:) + cl(idx_b1+4,:);
% a2l(idx_b,:) = -bl(idx_b1+2,:) - bl(idx_b1+3,:) + bl(idx_b1+4,:);
% a3l(idx_b,:) = -al(idx_b1+2,:) - al(idx_b1+3,:) + al(idx_b1+4,:);
% 
% idx_c = endpoints(idx+1,1,1) > endpoints(idx+1,1,2); % HDVEs with obtuse angle at Vertex 2
% idx_c1 = idx(idx_c.*idx ~= 0); 
% a1l(idx_c,:) = -cl(idx_c1+2,:) + cl(idx_c1+3,:) + cl(idx_c1+4,:);
% a2l(idx_c,:) = -bl(idx_c1+2,:) + bl(idx_c1+3,:) + bl(idx_c1+4,:);
% a3l(idx_c,:) = -al(idx_c1+2,:) + al(idx_c1+3,:) + al(idx_c1+4,:);
% 
% idx_d = abs(endpoints(idx,1,2)) <= dbl_eps; % HDVEs with right angle at Vertex 1
% idx_d1 = idx(idx_d.*idx ~= 0);
% a1l(idx_d,:) = -cl(idx_d1+3,:) + cl(idx_d1+4,:);
% a2l(idx_d,:) = -bl(idx_d1+3,:) + bl(idx_d1+4,:);
% a3l(idx_d,:) = -al(idx_d1+3,:) + al(idx_d1+4,:);
% 
% idx_e = abs(endpoints(idx+1,1,2) - endpoints(idx+1,1,1)) <= dbl_eps; % HDVEs with right angle at Vertex 3
% idx_e1 = idx(idx_e.*idx ~= 0);
% a1l(idx_e,:) = -cl(idx_e1+2,:) + cl(idx_e1+4,:);
% a2l(idx_e,:) = -bl(idx_e1+2,:) + bl(idx_e1+4,:);
% a3l(idx_e,:) = -al(idx_e1+2,:) + al(idx_e1+4,:);
% 
% a1l(dvetype == 2,:) = zeros(length(nonzeros(dvetype == 2)),3);
% a2l(dvetype == 2,:) = zeros(length(nonzeros(dvetype == 2)),3);
% a3l(dvetype == 2,:) = zeros(length(nonzeros(dvetype == 2)),3);
%% Transforming to global coordinates

v1 = [a1l; a2l; a3l; b1l; b2l; b3l];

dvenum = repmat(dvenum, 6, 1);

v2 = fcnSTARGLOB(v1, matROTANG(dvenum,1), matROTANG(dvenum,2), matROTANG(dvenum,3));

% v2 = v1;

a1 = v2(1:len,:);
a2 = v2(len+1:2*len,:);
a3 = v2(2*len+1:3*len,:);
b1 = v2(3*len+1:4*len,:);
b2 = v2(4*len+1:5*len,:);
b3 = v2(5*len+1:6*len,:);




