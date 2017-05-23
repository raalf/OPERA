function [PLEX, DVECT, ROTANG] = fcnTRITOLEXFS(P, DNORM)

% This function converts a 3-by-3-by-n matrix from global (X,Y,Z)
% to local (eta,xi,zeta) coordinates, where the depth n is any number of points
% and columns are (X,Y,Z), rows are vertex numbers 1,2,3
% Vertex number matters, as the first one is used as the origin and
% the third one is placed on the eta-axis
% T.D.K 2016-07-13. 212-230 KING ST E, TORONTO, ONTARIO, CANADA, M5A-1K5

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
DVECT(:,:,1) = global_x - ((dot(global_x,DVECT(:,:,3),2))./(sqrt(sum(DVECT(:,:,3).^2,2)))).*DVECT(:,:,3);

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
% UNFINISHED
[~,idx] = max(PLEX(:,1,:)); % Largest X value in the local coordinates



end