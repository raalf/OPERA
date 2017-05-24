function [PLEX, DVECT, ROTANG, matVSCOMB] = fcnTRITOLEXFS(P, DNORM)

% This function converts a nx3x3 matrix from global (X,Y,Z)
% to local (eta,xi,zeta) coordinates, where the depth n is vertex number
% and columns are (X,Y,Z), rows are HDVE number
% Vertex number matters, as the first one is used as the origin and
% the third one is placed on the eta-axis
% T.D.K 2016-07-13. 212-230 KING ST E, TORONTO, ONTARIO, CANADA, M5A-1K5

% matVSCOMB uses 1, 0, -1 to tell whether a vortex
% sheet is added or subtracted in the induction function
% which combines the vortex sheets on an HDVE
% 1 means the sheet is added, -1 means subtracted, 
% 0 means its perpendicular to the flow
% matVSCOMB(:,:,1) is for sheets going to positive eta infinity
% matVSCOMB(:,:,2) is for sheets going to positive xsi infinity

% Point 1 is the origin in local eta-xi reference frame, [0 0]

% Pre-allocating memory for a turbo-boost in performance
sp = size(P);
P2 = zeros(sp);

try
    PLEX = zeros(3,3,sp(3));
catch
    PLEX = zeros(3,3);
end
%zeta is always zero in the local plane of the hdve

%% Getting Roll/Pitch/Yaw
DVECT = zeros(sp(3),3,3);

% Normal
DVECT(:,:,3) = DNORM;

% Element "xsi" is aligned with global X (Change for rotors?)
global_x = repmat([1 0 0],sp(3),1);
DVECT(:,:,1) = global_x - repmat((dot(global_x,DVECT(:,:,3),2))./(sqrt(sum(DVECT(:,:,3).^2,2))),1,3).*DVECT(:,:,3);

% Element "eta" 
DVECT(:,:,2) = cross(DVECT(:,:,3),DVECT(:,:,1),2);

% Roll, pitch, yaw angles for global/local transformations
ROLL = -atan2(DVECT(:,2,3), DVECT(:,3,3));
PITCH = asin(DVECT(:,1,3));
YAW = atan2(DVECT(:,2,1), DVECT(:,1,1));

ROTANG(:,1) = ROLL;
ROTANG(:,2) = PITCH;
ROTANG(:,3) = YAW;

%% Local coordinates
temp_points = reshape(permute(P,[1 3 2]),[],size(P,2),1);
temp_points2 = fcnGLOBSTAR(temp_points, reshape(repmat(ROLL,1,3)',[],1),reshape(repmat(PITCH,1,3)',[],1), reshape(repmat(YAW,1,3)',[],1));
PLEX = permute(reshape(temp_points2',3,3,sp(3)),[2 1 3]);

%% Translating local coordinates towards origin
% Translating so that vertex 1 is at local origin (0,0)
PLEX = PLEX - repmat(PLEX(1,:,:),3,1,1);

%% Finding local leading and trailing edges of HDVEs (matVSCOMB)
a = (PLEX(2,:,:))./sqrt(sum(PLEX(2,:,:).^2,2));
c = (PLEX(3,:,:) - PLEX(2,:,:))./sqrt(sum((PLEX(3,:,:) - PLEX(2,:,:)).^2,2));
b = (PLEX(3,:,:))./sqrt(sum(PLEX(3,:,:).^2,2));

angle1 = atan2(a(:,2,:),a(:,1,:));
angle2 = atan2(c(:,2,:),c(:,1,:));
angle3 = atan2(b(:,2,:),b(:,1,:));

if isempty(nonzeros(angle3 > angle1))
    disp('Error in local reference frame edge numbering')
end
    
angle1(angle1 < 0) = angle1(angle1 < 0) + (2*pi);
angle2(angle2 < 0) = angle2(angle2 < 0) + (2*pi);
angle3(angle3 < 0) = angle3(angle3 < 0) + (2*pi);

matVSCOMB = zeros(sp(3),3,2);

%% For sheets extending in the eta direction to positive eta infinity
% Edge 1 and 2 are 1, Edge 3 is -1
idx = angle1 <= pi & angle2 < angle1;
len = length(nonzeros(idx));
matVSCOMB(idx,:,1) = repmat([1 1 -1],len,1);

% Edge 1 is 1, Edges 2 and 3 are -1
idx = angle1 <= pi & angle2 > angle1;
len = length(nonzeros(idx));
matVSCOMB(idx,:,1) = repmat([1 -1 -1],len,1);

% Edge 2 is 1, Edge 1 and 3 are -1
idx = angle1 >= pi & angle3 <= pi & angle1 - angle3 < pi;
len = length(nonzeros(idx));
matVSCOMB(idx,:,1) = repmat([-1 1 -1],len,1);

% Edge 1 and 3 are 1, Edge 2 is -1
idx = angle1 >= pi & angle3 <= pi & angle1 - angle3 > pi;
len = length(nonzeros(idx));
matVSCOMB(idx,:,1) = repmat([1 -1 1],len,1);

% Edge 2 and 3 are 1, Edge 1 is -1
idx = angle1 >= pi & angle2 <= pi & angle3 >= pi;
len = length(nonzeros(idx));
matVSCOMB(idx,:,1) = repmat([-1 1 1],len,1);

% Edge 3 is 1, Edge 1 and 2 are -1
idx = angle1 >= pi & angle2 >= pi & angle3 >= pi;
len = length(nonzeros(idx));
matVSCOMB(idx,:,1) = repmat([-1 -1 1],len,1);

% Special cases for above, when perpendicular to vortex sheet direction
idx = angle1 == 0 | angle1 == pi | angle1 == 2*pi;
len = length(nonzeros(idx));
matVSCOMB(idx,1,1) = zeros(len,1,1);

idx = angle2 == 0 | angle2 == pi | angle2 == 2*pi;
len = length(nonzeros(idx));
matVSCOMB(idx,2,1) = zeros(len,1,1);

idx = angle3 == 0 | angle3 == pi | angle3 == 2*pi;
len = length(nonzeros(idx));
matVSCOMB(idx,3,1) = zeros(len,1,1);

%% For sheets extending in the xsi direction to positive xsi infinity
% Edge 1 and 2 are 1, Edge 3 is -1
idx = angle1 <= pi/2 & (angle2 >= 3*pi/2 | angle2 < angle1);
len = length(nonzeros(idx));
matVSCOMB(idx,:,2) = repmat([1 1 -1],len,1);

% Edge 1 and 3 are -1, Edge2 is 1
idx = angle1 >= pi/2 & angle3 <= pi/2 & angle1 - angle3 < pi;
len = length(nonzeros(idx));
matVSCOMB(idx,:,2) = repmat([-1 1 -1],len,1);

% Edge 1 and 3 are 1, Edge 2 is -1
idx = angle1 >= pi/2 & angle1 <= 3*pi/2 & angle3 <= pi/2 & angle1 - angle3 > pi;
len = length(nonzeros(idx));
matVSCOMB(idx,:,2) = repmat([1 -1 1],len,1);

% Edge 1 is -1, Edges 2 and 3 are 1
idx = angle1 >= pi/2 & angle3 >= pi/2 & angle1 <= 3*pi/2 & angle3 <= 3*pi/2 & angle2 <= pi/2;
len = length(nonzeros(idx));
matVSCOMB(idx,:,2) = repmat([-1 1 1],len,1);

% Edge 1 and 2 are -1, Edge 3 is 1
idx = angle1 >= pi/2 & angle3 >= pi/2 & angle1 <= 3*pi/2 & angle3 <= 3*pi/2 & angle2 >= pi/2;
len = length(nonzeros(idx));
matVSCOMB(idx,:,2) = repmat([-1 -1 1],len,1);

% Edge 1 and 3 are 1, Edge 2 is -1
idx = angle1 >= 3*pi/2 & angle3 <= 3*pi/2 & angle1 - angle3 < pi;
len = length(nonzeros(idx));
matVSCOMB(idx,:,2) = repmat([1 -1 1],len,1);

% Edge 1 and 3 are -1, Edge 2 is 1
idx = angle1 >= 3*pi/2 & angle3 <= 3*pi/2 & angle3 >= pi/2 & angle1 - angle3 > pi;
len = length(nonzeros(idx));
matVSCOMB(idx,:,2) = repmat([-1 1 -1],len,1);

% Edge 3 and 2 are 1, Edge 1 is -1
idx = angle1 >= 3*pi/2 & angle3 <= pi/2 & angle2 >= pi/2;
len = length(nonzeros(idx));
matVSCOMB(idx,:,2) = repmat([-1 1 1],len,1);

% Edge 1 and 2 are -1, Edge 3 is 1
idx = angle1 >= 3*pi/2 & angle3 <= pi/2 & angle2 <= pi/2;
len = length(nonzeros(idx));
matVSCOMB(idx,:,2) = repmat([-1 -1 1],len,1);

% Edge 1 is 1, Edges 2 and 3 are -1
idx = angle1 >= 3*pi/2 & angle3 >= 3*pi/2;
len = length(nonzeros(idx));
matVSCOMB(idx,:,2) = repmat([1 -1 -1],len,1);

% Special cases for above, when perpendicular to vortex sheet direction
idx = angle1 == pi/2 | angle1 == 3*pi/2;
len = length(nonzeros(idx));
matVSCOMB(idx,1,2) = zeros(len,1,1);

idx = angle2 == pi/2 | angle2 == 3*pi/2;
len = length(nonzeros(idx));
matVSCOMB(idx,2,2) = zeros(len,1,1);

idx = angle3 == pi/2 | angle3 == 3*pi/2;
len = length(nonzeros(idx));
matVSCOMB(idx,3,2) = zeros(len,1,1);

end