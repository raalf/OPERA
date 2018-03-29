function [xsi] = fcnSTARGLOB(points, rotation_angles)
% Transforms a point from local to global reference frame
% Input in RADIANS

% INPUT:
%   points - n x 3 matrix of (x,y,z) points in local coordinates
%   rotation_angles - n x 3 matrix, where colums are (roll, pitch, yaw)
% OUTPUT:
%   xsi - n x 3 matrix of "points" in global reference frame

len = length(points(:,1));

cnu = cos(rotation_angles(:,1));
snu = sin(rotation_angles(:,1));
ceps = cos(rotation_angles(:,2));
seps = sin(rotation_angles(:,2));
cpsi = cos(rotation_angles(:,3));
spsi = sin(rotation_angles(:,3));

xsi = zeros(len,3);
xsi(:,1) = points(:,1).*(cpsi.*ceps) + points(:,2).*(-spsi.*ceps) + points(:,3).*(seps);
xsi(:,2) = points(:,1).*(cpsi.*seps.*snu+spsi.*cnu) + points(:,2).*(-spsi.*seps.*snu+cpsi.*cnu) + points(:,3).*(-ceps.*snu);
xsi(:,3) = points(:,1).*(-cpsi.*seps.*cnu+spsi.*snu) + points(:,2).*(spsi.*seps.*cnu+cpsi.*snu) + points(:,3).*(ceps.*cnu);

end

