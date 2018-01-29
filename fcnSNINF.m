function [infl] = fcnSNINF(qm, qn, c3, z, h, N, flagPLOT, ii, pb, ax2, loc_rot)

c2 = (qn - qm);
c2 = c2./sqrt(sum(c2.^2,2));
c1 = cross(c2,c3);

% idx_flip = sign(dot(qm,c1, 2)) < 0 & sign(dot(qm,c2,2)) < 0;
idx_flip = sign(dot(qm,c1, 2)) < 0;
c2(idx_flip,:) = c2(idx_flip,:).*-1;
c1(idx_flip,:) = c1(idx_flip,:).*-1;
temp = qm;
qm(idx_flip,:) = qn(idx_flip,:);
qn(idx_flip,:) = temp(idx_flip,:);


ori = dot(c3, cross(qm, qn), 2)./abs(dot(c3, cross(qm, qn), 2));
ori(isnan(ori)) = 1;

comb = dot(N, qm, 2);
comb = comb./abs(comb);
comb(isnan(comb)) = 0;

% len = length(h);
% loc_loc_rot = [-atan2(c3(:,2), c3(:,3)), asin(c3(:,1)), atan2(c1(:,2), c1(:,1))];
% 
% qm = fcnGLOBSTAR(qm, loc_loc_rot(:,1), loc_loc_rot(:,2), loc_loc_rot(:,3));
% qn = fcnGLOBSTAR(qn, loc_loc_rot(:,1), loc_loc_rot(:,2), loc_loc_rot(:,3));
% c1 = fcnGLOBSTAR(c1, loc_loc_rot(:,1), loc_loc_rot(:,2), loc_loc_rot(:,3));
% c2 = fcnGLOBSTAR(c2, loc_loc_rot(:,1), loc_loc_rot(:,2), loc_loc_rot(:,3));
% c3 = fcnGLOBSTAR(c3, loc_loc_rot(:,1), loc_loc_rot(:,2), loc_loc_rot(:,3));

a = dot(qm, c1, 2);
l1 = dot(qm, c2, 2);
l2 = dot(qn, c2, 2);

phi = [];
phi(:,1) = atan2d(l1, a);
phi(:,2) = atan2d(l2, a);

if flagPLOT == 1 && size(h,1) == 1
   temp_plt_s2;    
end

[H00, H10, H01, H20, H11, H02] = fcnHINTEGRAL(a,h,l1,l2);

a1c = [-z.*H01.*c2(:,1), -z.*H01.*c2(:,2), -H02.*c3(:,3)];
a2c = [-z.*H00.*c2(:,1), -z.*H00.*c2(:,2), -H01.*c3(:,3)];
b1c = [z.*H10.*c1(:,1), z.*H10.*c1(:,2), H20.*c3(:,3)];
b2c = [z.*H00.*c1(:,1), z.*H00.*c1(:,2), H10.*c3(:,3)];
c2c = [z.*H01.*c1(:,1) - z.*H10.*c2(:,1), z.*H01.*c1(:,2) - z.*H10.*c2(:,2), c3(:,3).*0];
c3c = c2c.*0;

% % Orientation of the triangle (determines whether to add or subtract influence)
% % First, set all orientations the same, then apply comb to add or subtract

% Order: A1 A2 B1 B2 C2 C3
infl = [reshape(a1c',3,1,[]) reshape(a2c',3,1,[]) reshape(b1c',3,1,[]) reshape(b2c',3,1,[]) reshape(c2c',3,1,[]) reshape(c3c',3,1,[])];

infl = infl.*repmat(reshape(ori,1,1,[]),3,6,1).*repmat(reshape(comb,1,1,[]),3,6,1);
% infl = infl.*repmat(reshape(ori,1,1,[]),3,6,1);

% dvenum = reshape(repmat([1:len]',1,6,1)',[],1,1);
% % Rotate from local c to local b
% infl = fcnSTARGLOB(reshape(permute(infl,[2 3 1]),[],3,1), loc_loc_rot(dvenum,1), loc_loc_rot(dvenum,2), loc_loc_rot(dvenum,3));
% infl = reshape(infl',3,6,[]);
end