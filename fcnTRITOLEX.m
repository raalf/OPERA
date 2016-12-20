function [PLEX, DVECT, ROTANG] = fcnTRITOLEX(P, DNORM)

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

% Below we create vectors from Point 1 to Points 2 and 3 in local x-y-z
P2(2,:,:) = P(2,:,:) - P(1,:,:);
P2(3,:,:) = P(3,:,:) - P(1,:,:);

% Lengths of all three sides
t = sqrt(sum(abs(P2).^2,2));
t1 = t(2:3:end,:,:); % Side 1
t3 = t(3:3:end,:,:); % Side 3

DVECT(:,:,1) = permute(P2(3,:,:)./repmat(t3,1,sp(2),1), [3 2 1]);
DVECT(:,:,2) = cross(DVECT(:,:,1), -DNORM); % DNORM is negative to ensure we get the proper 
DVECT(:,:,3) = DNORM;
%axang2rotm

% Side 2 is different, because it is between points 2 and 3, and that vector is not yet defined
t2 = sqrt(sum(abs(P2(3:3:end,:,:)-P2(2:3:end,:,:)).^2,2));

% Magnitude from Point 1 to Point 3 is the eta component of the local eta-xi reference frame
% (xi is zero because Point 3 lies on the eta-axis)
PLEX(3,1,:) = t3;

% Second point calculation using cosine law
% This seems to work fine when gam > 90 degrees, I'm not sure if this could be a problem later on
gam = acos( (t1.^2 + t3.^2 - t2.^2)./(2.*t1.*t3));
PLEX(2,1:2,:) = [t1.*cos(gam) t1.*sin(gam)];

ROLL = -atan2(DVECT(:,2,3), DVECT(:,3,3));
PITCH = asin(DVECT(:,1,3));
YAW = acos(dot(DVECT(:,:,2),repmat([0 1 0],sp(2),1,1),2));

ROTANG(:,1) = ROLL;
ROTANG(:,2) = PITCH;
ROTANG(:,3) = YAW;

end