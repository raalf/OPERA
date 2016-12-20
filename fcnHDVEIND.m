function [a1, a2, b1, b2, c3] = fcnHDVEIND(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX, matROTANG)

dvenum = reshape(dvenum, [], 1, 1); % Ensuring dvenum is a column vector

len = length(dvenum);
dbl_eps = 1e-14;

fp = fcnTOLOC(dvenum, fpg, matDVE, matDVECT, matVLST);

% fp = fcnGLOBSTAR(fpg - matVLST(matDVE(dvenum,1,1),:), matROTANG(:,1), matROTANG(:,2), matROTANG(:,3));



endpoints = zeros(len*5, 3, 2); % Five function calls per DVE
phi = zeros(len*5,1);
yaw = zeros(len*5,1);

k = zeros(len*5,1)+1; % Until this is figured out for HDVE, k = 1

idx = [1:5:len*5]';

%% Finding velocities induced by each vortex sheet
% Always input left point, right point of vortex sheet LE

% Left to right
% First edge
endpoints(idx,:,1:2) = reshape(permute(matPLEX(1:2,:,dvenum),[3 2 1]), [],3,2);
yaw(idx) = pi/2;
phi(idx) = abs(atan(endpoints(idx,1,2)./endpoints(idx,2,2)));

% Second edge
endpoints(idx+1,:,1:2) = reshape(permute(matPLEX(3:-1:2,:,dvenum),[3 2 1]), [],3,2);
yaw(idx+1) = pi/2;
phi(idx+1) = abs(atan((endpoints(idx+1,1,1)-endpoints(idx+1,1,2))./endpoints(idx+1,2,2)));

% Up to down
% First edge
endpoints(idx+2,:,1:2) = reshape(permute(matPLEX(1:2,:,dvenum),[3 2 1]), [],3,2);
yaw(idx+2) = 0;
phi(idx+2) = abs(atan(endpoints(idx+2,2,2)./endpoints(idx+2,1,2)));

% Second edge
endpoints(idx+3,:,1:2) = reshape(permute(matPLEX(2:3,:,dvenum),[3 2 1]), [],3,2);
yaw(idx+3) = 0;
phi(idx+3) = abs(atan(endpoints(idx+3,2,1)./(endpoints(idx+3,1,2)-endpoints(idx+3,1,1))));

% Third edge
endpoints(idx+4,:,1:2) = reshape(permute(matPLEX(1:2:3,:,dvenum),[3 2 1]), [],3,2);
yaw(idx+4) = 0;
phi(idx+4) = 0; % Leading edge of Edge 3 is always on eta-axis

%% Case with obtuse angle at Vertex 1
idx_a = endpoints(idx+2,1,2) < -dbl_eps; % HDVEs with obtuse angle at Vertex 1 (eta < 0)
idx_a1 = idx(idx_a.*idx ~= 0); % Index of first edge of HDVEs for obtuse angle at Vertex 1

if ~isempty(idx_a1) == 1
    % In this special case, Edge 2 is the longest edge (idx+3). What we need to do is extend the spans
    % of the 1st and 3rd vortex sheet to match the span of the second. Then we need to evaluate HVSIND only
    % on the range where the edge actually exists
    
    % Swapping endpoints so this vortex sheet goes from left to right (Edge 1)
    temp1 = endpoints(idx_a1+2,:,1);
    endpoints(idx_a1+2,:,1) = endpoints(idx_a1+2,:,2);
    endpoints(idx_a1+2,:,2) = temp1;
    
    % Adjusting phi so it is the correct sign
    phi(idx_a1) = -phi(idx_a1);
    phi(idx_a1+1) = -phi(idx_a1+1);
    
end


%% Case with acute angles at Vertices 1 and 3
idx_b = endpoints(idx+2,1,2) > dbl_eps & endpoints(idx+2,1,2) < (endpoints(idx+4,1,2) - dbl_eps); % HDVEs with acute (or right) angle at Vertex 1 & Vertex 3
idx_b1 = idx(idx_b.*idx ~= 0);

if ~isempty(idx_b1) == 1
    % All three vortex sheets should be in the same direction already (left to right)
    
    % Adjusting phi so it is the correct sign
    phi(idx_b1+1) = -phi(idx_b1+1);
    phi(idx_b1+2) = -phi(idx_b1+2);
    
end
%% Case with obtuse angle at Vertex 3
idx_c = endpoints(idx+2,1,2) > (endpoints(idx+4,1,2) + dbl_eps); % HDVEs with obtuse angle at Vertex 2
idx_c1 = idx(idx_c.*idx ~= 0);

if ~isempty(idx_c1) == 1
    % Swapping endpoints so this vortex sheet goes from left to right (Edge 2)
    temp1 = endpoints(idx_c1+3,:,1);
    endpoints(idx_c1+3,:,1) = endpoints(idx_c1+3,:,2);
    endpoints(idx_c1+3,:,2) = temp1;

    % Adjusting phi so it is the correct sign
    phi(idx_c1+2) = -phi(idx_c1+2);
    phi(idx_c1+3) = -phi(idx_c1+3);
   
end

%% Right angle at Vertex 1
idx_d = abs(endpoints(idx,1,2)) <= dbl_eps; % HDVEs with right angle at Vertex 2
idx_d1 = idx(idx_d.*idx ~= 0);

if ~isempty(idx_d1) == 1
    % Adjusting phi so it is the correct sign
    phi(idx_d1) = 0;
    phi(idx_d1+1) = -phi(idx_d1+1);
end
%% Right angle at Vertex 3
idx_e = abs(endpoints(idx+1,1,2) - endpoints(idx+1,1,1)) <= dbl_eps;
idx_e1 = idx(idx_e.*idx ~= 0);

if ~isempty(idx_e1) == 1
    % Adjusting phi so it is the correct sign
    phi(idx_e1+1) = 0;
    phi(idx_e1+2) = -phi(idx_e1+2);
end
%% Running VSIND

phi = -phi;

fpl = reshape(repmat(fp,1,5,1)',3,[],1)';

[al, bl, cl] = fcnVSIND(endpoints, phi, yaw, fpl, k);

a1l = zeros(len,3);
a2l = zeros(len,3);
a3l = zeros(len,3);

b1l = zeros(len,3);
b2l = zeros(len,3);
b3l = zeros(len,3);

%% Summing the velocities of all five sheets

% Subtracting second vortex sheet from first in the left-to-right direction
b1l = cl(idx,:) - cl(idx+1,:);
b2l = bl(idx,:) - bl(idx+1,:);
b3l = al(idx,:) - al(idx+1,:);

% Subtracting the three vortex sheets in the up-to-down direction based on the shape of the triangle

idx_a = endpoints(idx,1,2) < 0; % HDVEs with obtuse angle at Vertex 1 (eta < 0)
idx_a1 = idx(idx_a.*idx ~= 0); % Index of first edge of HDVEs for obtuse angle at Vertex 1
a1l(idx_a,:) = cl(idx_a1+2,:) - cl(idx_a1+3,:) + cl(idx_a1+4,:);
a2l(idx_a,:) = bl(idx_a1+2,:) - bl(idx_a1+3,:) + bl(idx_a1+4,:);
a3l(idx_a,:) = al(idx_a1+2,:) - al(idx_a1+3,:) + al(idx_a1+4,:);

idx_b = endpoints(idx,1,2) >= 0 & endpoints(idx+1,1,2) <= endpoints(idx+1,1,1); % HDVEs with acute (or right) angle at Vertex 1 & Vertex 2
idx_b1 = idx(idx_b.*idx ~= 0);
a1l(idx_b,:) = -cl(idx_b1+2,:) - cl(idx_b1+3,:) + cl(idx_b1+4,:);
a2l(idx_b,:) = -bl(idx_b1+2,:) - bl(idx_b1+3,:) + bl(idx_b1+4,:);
a3l(idx_b,:) = -al(idx_b1+2,:) - al(idx_b1+3,:) + al(idx_b1+4,:);

idx_c = endpoints(idx+1,1,2) > endpoints(idx+1,1,1); % HDVEs with obtuse angle at Vertex 2
idx_c1 = idx(idx_c.*idx ~= 0); 
a1l(idx_c,:) = -cl(idx_c1+2,:) + cl(idx_c1+3,:) + cl(idx_c1+4,:);
a2l(idx_c,:) = -bl(idx_c1+2,:) + bl(idx_c1+3,:) + bl(idx_c1+4,:);
a3l(idx_c,:) = -al(idx_c1+2,:) + al(idx_c1+3,:) + al(idx_c1+4,:);

idx_d = abs(endpoints(idx,1,2)) <= dbl_eps; % HDVEs with right angle at Vertex 1
idx_d1 = idx(idx_d.*idx ~= 0);
a1l(idx_d,:) = -cl(idx_d1+3,:) + cl(idx_d1+4,:);
a2l(idx_d,:) = -bl(idx_d1+3,:) + bl(idx_d1+4,:);
a3l(idx_d,:) = -al(idx_d1+3,:) + al(idx_d1+4,:);

idx_e = abs(endpoints(idx+1,1,2) - endpoints(idx+1,1,1)) <= dbl_eps; % HDVEs with right angle at Vertex 2
idx_e1 = idx(idx_e.*idx ~= 0);
a1l(idx_e,:) = -cl(idx_e1+2,:) + cl(idx_e1+4,:);
a2l(idx_e,:) = -bl(idx_e1+2,:) + bl(idx_e1+4,:);
a3l(idx_e,:) = -al(idx_e1+2,:) + al(idx_e1+4,:);


%% Transforming to global coordinates

v1 = [a1l; a2l; b1l; b2l; a3l+b3l];

% v2 = fcnSTARGLOB(v1, matROTANG(:,1), matROTANG(:,2), matROTANG(:,3));
v2 = fcnTOGLOB(reshape(repmat(dvenum,1,5)',[],1), v1, matDVE, matDVECT, matVLST);
% v2 = fcnROTVECT(repmat(dvenum,5,1,1), v1, matDVECT);
% v2 = v1

a1 = v2(1:len,:);
a2 = v2(len+1:2*len,:);
b1 = v2(2*len+1:3*len,:);
b2 = v2(3*len+1:4*len,:);
c3 = v2(4*len+1:5*len,:);



