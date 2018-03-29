function [PLEX, DVECT, ROTANG] = fcnTRITOLEX(P, DNORM, matCENTER)
% This function converts a nx3x3 matrix from global (X,Y,Z)
% to local (eta,xi,zeta) coordinates, where the depth n is vertex number
% and columns are (X,Y,Z), rows are HDVE number
% Vertex number matters maybe?
% The local origin is set to the incenter of the triangle
% T.D.K 2016-07-13. 230 KING ST E, TORONTO, ONTARIO, CANADA, M5A-1K5

% matVSCOMB uses 1, 0, -1 to tell whether a vortex
% sheet is added or subtracted in the induction function
% which combines the vortex sheets on an HDVE
% 1 means the sheet is added, -1 means subtracted, 
% 0 means its perpendicular to the flow
% matVSCOMB(:,:,1) is for sheets going to positive eta infinity
% matVSCOMB(:,:,2) is for sheets going to positive xsi infinity
% MODIFIED
% T.D.K 2017-06-05, 907 E 13TH AVE, DENVER, COLORADO, USA, 80218

% Gutted, rewrote, origin is vertex 1 and one side is aligned with local
% eta axis
% T.D.K 2018-03-17. 230 KING ST E, TORONTO, ONTARIO, CANADA, M5A-1K5

% Pre-allocating memory for a turbo-boost in performance
sp = size(P);
P2 = zeros(sp);

try
    PLEX = zeros(3,3,sp(3));
catch
    PLEX = zeros(3,3);
end

%% Getting Roll/Pitch/Yaw
if size(sp,2) == 2
    sp = [sp 1];
end

DVECT = zeros(sp(3),3,3);

% Zeta (Normal)
DVECT(:,:,3) = DNORM;

% Eta (think of it as local y-axis)
DVECT(:,:,2) = permute((P(2,:,:) - P(1,:,:))./sqrt(sum((P(2,:,:) - P(1,:,:)).^2,2)), [3 2 1]);

% Xsi
DVECT(:,:,1) = cross(DVECT(:,:,2),DVECT(:,:,3),2);

% Roll, pitch, yaw angles for global/local transformations
% [YAW, PITCH, ROLL] = dcm2angle(permute(DVECT, [3 2 1]), 'ZXY');
[ROLL, PITCH, YAW] = dcm2angle(permute(DVECT, [3 2 1]), 'XYZ', 'ZeroR3');

ROTANG(:,1) = ROLL;
ROTANG(:,2) = PITCH;
ROTANG(:,3) = YAW;

%% Local coordinates
temp_points = reshape(permute(P,[1 3 2]),[],size(P,2),1);
temp_points2 = fcnGLOBSTAR(temp_points, reshape(repmat(ROLL,1,3)',[],1),reshape(repmat(PITCH,1,3)',[],1), reshape(repmat(YAW,1,3)',[],1));
PLEX = permute(reshape(temp_points2',3,3,sp(3)),[2 1 3]);

%% Translating local coordinates towards origin
% Translating so that local origin is at the triangle incenter (matCENTER)
temp_center = fcnGLOBSTAR(matCENTER, ROLL, PITCH, YAW);
PLEX = PLEX - repmat(reshape(temp_center',1,3,[]),3,1,1);

%% Finding local leading and trailing edges of HDVEs (matVSCOMB)
matVSCOMB = zeros(sp(3),3,2);

eta = zeros(1,3,sp(3));
xsi = zeros(1,3,sp(3));

% Edge vectors by going around in order
edge = [PLEX(1,:,:) - PLEX(2,:,:); PLEX(2,:,:) - PLEX(3,:,:); PLEX(3,:,:) - PLEX(1,:,:)];

% Finding edge normals and normalizing
edge_normal = [edge(:,2,:) -edge(:,1,:) edge(:,3,:)];
edge_normal = edge_normal./(sqrt(sum(edge_normal(:,1,:).^2 + edge_normal(:,2,:).^2,2)));

check_normal = dot(edge_normal(1,:,:), edge(3,:,:));
check_normal(check_normal > 0) = 1; % Here we need to flip the normals
check_normal(check_normal <= 0) = 0;

edge_normal(:,:,check_normal == 1) = -edge_normal(:,:,check_normal ==1);

eta(permute(edge_normal(:,1,:) > 0, [2 1 3])) = -1;
eta(permute(edge_normal(:,1,:) < 0, [2 1 3])) = 1;
matVSCOMB(:,:,1) = reshape(permute(eta, [2 1 3]), size(eta, 2), [])';

xsi(permute(edge_normal(:,2,:) > 0, [2 1 3])) = -1;
xsi(permute(edge_normal(:,2,:) < 0, [2 1 3])) = 1;
matVSCOMB(:,:,2) = reshape(permute(xsi, [2 1 3]), size(xsi, 2), [])';


end