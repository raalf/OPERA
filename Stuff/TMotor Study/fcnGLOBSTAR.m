function [xsi] = fcnGLOBSTAR(points, rotation_angles)
% Transforms a point from global to local reference frame
% Input in RADIANS

% INPUT:
%   points - n x 3 matrix of (x,y,z) points in global coordinates
%   rotation_angles - n x 3 matrix, where colums are (roll, pitch, yaw)
% OUTPUT:
%   xsi - n x 3 matrix of "points" in local reference frame

len = length(points(:,1));

cnu = cos(rotation_angles(:,1));
snu = sin(rotation_angles(:,1));
ceps = cos(rotation_angles(:,2));
seps = sin(rotation_angles(:,2));
cpsi = cos(rotation_angles(:,3));
spsi = sin(rotation_angles(:,3));

xsi = zeros(len,3);
xsi(:,1) = points(:,1).*(cpsi.*ceps) + points(:,2).*(cpsi.*seps.*snu+spsi.*cnu) + points(:,3).*(-cpsi.*seps.*cnu+spsi.*snu);
xsi(:,2) = points(:,1).*(-spsi.*ceps) + points(:,2).*(-spsi.*seps.*snu+cpsi.*cnu) + points(:,3).*(spsi.*seps.*cnu+cpsi.*snu);
xsi(:,3) = points(:,1).*(seps) + points(:,2).*(-ceps.*snu) + points(:,3).*(ceps.*cnu);

end

