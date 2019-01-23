function [PLEX, DVECT, ROTANG] = fcnTRITOLEX(P, DNORM, matCONTROL, TYPE, matWETA)
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

if strcmpi(TYPE, 'SURFACE')
    % Element "xsi" is aligned with global X (Change for rotors?)
    global_x = repmat([-1 0 0],sp(3),1);
    DVECT(:,:,2) = global_x - repmat((dot(global_x,DVECT(:,:,3),2))./(sqrt(sum(DVECT(:,:,3).^2,2))),1,3).*DVECT(:,:,3);
else
    idx = reshape(matWETA,1,1,[]);
    DVECT(matWETA == 1,:,2) = reshape(P(2,:,idx == 1) - P(1,:,idx == 1),3,[],1)';
    DVECT(matWETA == 2,:,2) = reshape(P(2,:,idx == 2) - P(3,:,idx == 2),3,[],1)';
end

% Element "eta" 
DVECT(:,:,1) = cross(DVECT(:,:,2),DVECT(:,:,3),2);

% % Eta (think of it as local y-axis)
% DVECT(:,:,2) = permute((P(2,:,:) - P(1,:,:))./sqrt(sum((P(2,:,:) - P(1,:,:)).^2,2)), [3 2 1]);
% 
% % Xsi
% DVECT(:,:,1) = cross(DVECT(:,:,2),DVECT(:,:,3),2);

% Roll, pitch, yaw angles for global/local transformations
% [YAW, PITCH, ROLL] = dcm2angle(permute(DVECT, [3 2 1]), 'ZXY');
[ROLL, PITCH, YAW] = dcm2angle(permute(DVECT, [3 2 1]), 'XYZ', 'ZeroR3');

ROTANG(:,1) = ROLL;
ROTANG(:,2) = PITCH;
ROTANG(:,3) = YAW;

%% Local coordinates
temp_points = reshape(permute(P,[1 3 2]),[],size(P,2),1);
temp_points2 = fcnGLOBSTAR(temp_points, [reshape(repmat(ROLL,1,3)',[],1) reshape(repmat(PITCH,1,3)',[],1) reshape(repmat(YAW,1,3)',[],1)]);
PLEX = permute(reshape(temp_points2',3,3,sp(3)),[2 1 3]);

%% Translating local coordinates towards origin
% Translating so that local origin is at the triangle incenter (matCENTER)
temp_center = fcnGLOBSTAR(matCONTROL, ROTANG);
PLEX = PLEX - repmat(reshape(temp_center',1,3,[]),3,1,1);

end