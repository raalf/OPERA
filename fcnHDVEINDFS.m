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
yaw(idx) = 0;
phi(idx) = atan((endpoints(idx,1,2)-endpoints(idx,1,1))./(endpoints(idx,2,2)-endpoints(idx,2,1)));

% Second edge
endpoints(idx+1,:,1:2) = reshape(permute(matPLEX(2:3,:,dvenum),[3 2 1]), [],3,2);
yaw(idx+1) = 0;
phi(idx+1) = atan((endpoints(idx+1,1,2)-endpoints(idx+1,1,1))./(endpoints(idx+1,2,2)-endpoints(idx+1,2,1)));

% Third edge
endpoints(idx+2,:,1:2) = reshape(permute(matPLEX(1:2:3,:,dvenum),[3 2 1]), [],3,2);
yaw(idx+2) = 0;
phi(idx+2) = atan((endpoints(idx+2,1,2)-endpoints(idx+2,1,1))./(endpoints(idx+2,2,2)-endpoints(idx+2,2,1)));

% Up to down. Sheet LE is on the bottom, sheet extends to infinity upwards.
% First edge
endpoints(idx+3,:,1:2) = reshape(permute(matPLEX(1:2,:,dvenum),[3 2 1]), [],3,2);
yaw(idx+3) = pi/2;
phi(idx+3) = -atan((endpoints(idx+3,2,2)-endpoints(idx+3,2,1))./(endpoints(idx+3,1,2)-endpoints(idx+3,1,1)));

% Second edge
endpoints(idx+4,:,1:2) = reshape(permute(matPLEX(2:3,:,dvenum),[3 2 1]), [],3,2);
yaw(idx+4) = pi/2;
phi(idx+4) = -atan((endpoints(idx+4,2,2)-endpoints(idx+4,2,1))./(endpoints(idx+4,1,2)-endpoints(idx+4,1,1)));

% Third edge
endpoints(idx+5,:,1:2) = reshape(permute(matPLEX(1:2:3,:,dvenum),[3 2 1]), [],3,2);
yaw(idx+5) = pi/2;
phi(idx+5) = -atan((endpoints(idx+5,2,2)-endpoints(idx+5,2,1))./(endpoints(idx+5,1,2)-endpoints(idx+5,1,1))); % Leading edge of Edge 3 is always on eta-axis

%% Running VSIND

fpl = reshape(repmat(fp,1,6,1)',3,[],1)';

[al, bl, cl] = fcnVSIND(endpoints, phi, yaw, fpl, k); 

idx_null = abs(phi) == pi/2;

bl(idx_null,:) = zeros(length(nonzeros(idx_null)),3);
cl(idx_null,:) = zeros(length(nonzeros(idx_null)),3);

a1l = zeros(len,3);
a2l = zeros(len,3);
a3l = zeros(len,3);

%% Summing the velocities of all five sheets

% RIGHT-TO-LEFT
b1l = repmat(matVSCOMB(dvenum,1,1),1,3).*cl(idx,:) + repmat(matVSCOMB(dvenum,2,1),1,3).*cl(idx+1,:) + repmat(matVSCOMB(dvenum,3,1),1,3).*cl(idx+2,:);
b2l = repmat(matVSCOMB(dvenum,1,1),1,3).*bl(idx,:) + repmat(matVSCOMB(dvenum,2,1),1,3).*bl(idx+1,:) + repmat(matVSCOMB(dvenum,3,1),1,3).*bl(idx+2,:);
b3l = repmat(matVSCOMB(dvenum,1,1),1,3).*al(idx,:) + repmat(matVSCOMB(dvenum,2,1),1,3).*al(idx+1,:) + repmat(matVSCOMB(dvenum,3,1),1,3).*al(idx+2,:);

% DOWN-TO-UP (negative because matVSCOMB assumes sheet goes left to right, to positive eta infinity)
% matVSCOMB(:,:,2) = matVSCOMB(:,:,2)*-1;
a1l = repmat(matVSCOMB(dvenum,1,2),1,3).*cl(idx+3,:) + repmat(matVSCOMB(dvenum,2,2),1,3).*cl(idx+4,:) + repmat(matVSCOMB(dvenum,3,2),1,3).*cl(idx+5,:);
a2l = repmat(matVSCOMB(dvenum,1,2),1,3).*bl(idx+3,:) + repmat(matVSCOMB(dvenum,2,2),1,3).*bl(idx+4,:) + repmat(matVSCOMB(dvenum,3,2),1,3).*bl(idx+5,:);
a3l = repmat(matVSCOMB(dvenum,1,2),1,3).*al(idx+3,:) + repmat(matVSCOMB(dvenum,2,2),1,3).*al(idx+4,:) + repmat(matVSCOMB(dvenum,3,2),1,3).*al(idx+5,:);

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




