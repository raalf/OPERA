function v2 = fcnROTVECT(dvenum, v1, DVECT)
% Takes a vector in local DVE reference frame and rotates it 
% so it is aligned with the global X-Y-Z axis system.

% INPUT:
%   dvenum - n x 1 vector of DVE number to rotate with respect to
%   v1 - n x 3 matrix of points to be rotated
%   DVECT

% OUTPUT:
%   v2 - v1 rotated so it is aligned with global axis

% T.D.K 2016-09-08 745-55 RIVER OAKS PLACE, SAN JOSE, CALIFORNIA, USA 95134

len = length(dvenum);

alpha_r = acos(dot(DVECT(dvenum,:,1),repmat([1 0 0],len,1,1),2));
beta_r = acos(dot(DVECT(dvenum,:,2),repmat([0 1 0],len,1,1),2));
gamma_r = acos(dot(DVECT(dvenum,:,3),repmat([0 0 1],len,1,1),2));

R_alpha = repmat([0 0 0],3,1,len);
R_alpha(1,1,:) = ones(len,1);
R_alpha(2,2,:) = cos(alpha_r);
R_alpha(2,3,:) = sin(alpha_r);
R_alpha(3,2,:) = -sin(alpha_r);
R_alpha(3,3,:) = cos(alpha_r);

R_beta = repmat([0 0 0],3,1,len);
R_beta(1,1,:) = cos(beta_r);
R_beta(1,3,:) = -sin(beta_r);
R_beta(2,2,:) = ones(len,1);
R_beta(3,1,:) = sin(beta_r);
R_beta(3,3,:) = cos(beta_r);

R_gamma = repmat([0 0 0],3,1,len);
R_gamma(1,1,:) = cos(gamma_r);
R_gamma(1,2,:) = sin(gamma_r);
R_gamma(2,1,:) = -sin(gamma_r);
R_gamma(2,2,:) = cos(gamma_r);
R_gamma(3,3,:) = ones(len,1);

v1t = repmat(permute(v1,[3 2 1]),3,1,1);

Rot = R_alpha.*R_beta.*R_gamma;

v2 = permute(sum(v1t.*Rot,2),[2 1 3]);
v2 = reshape(permute(v2,[3 1 2]),[],3,1);

end