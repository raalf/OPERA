function [infl] = fcnSNINF(qm, qn, c3, b, z, h)

c2 = (qn - qm);
c2 = c2./sqrt(sum(c2.^2,2));
c1 = cross(c2,c3);

% % WHAAAAAAAAAAAAAAAAAAAAAAAAT
% c1 = b(:,:,1);
% c2 = b(:,:,2);

a = dot(qm, c1, 2);
l1 = dot(qm, c2, 2);
l2 = dot(qn, c2, 2);

[H00, H10, H01, H20, H11, H02] = fcnHINTEGRAL_2(a,h,l1,l2);

a1 = [-z.*H01.*c2(:,1), -z.*H01.*c2(:,2), -H02.*c3(:,3)];
a2 = [-z.*H00.*c2(:,1), -z.*H00.*c2(:,2), -H01.*c3(:,3)];

b1 = [z.*H10.*c1(:,1), z.*H10.*c1(:,2), H20.*c3(:,3)];
b2 = [z.*H00.*c1(:,1), z.*H00.*c1(:,2), H10.*c3(:,3)];

c2 = [z.*H01.*c1(:,1) - z.*H10.*c2(:,1), z.*H01.*c1(:,2) - z.*H10.*c2(:,2), H11.*c3(:,3) - H11.*c3(:,3)];

% Orientation of the triangle (determines whether to add or subtract influence)
% First, set all orientations the same, then apply comb to add or subtract
ori = dot(c3, cross(qm, qn), 2)./abs(dot(c3, cross(qm, qn), 2));
ori(isnan(ori)) = 1;

% Order: A1 A2 B1 B2 C2
infl = [reshape(a1',3,1,[]) reshape(a2',3,1,[]) reshape(b1',3,1,[]) reshape(b2',3,1,[]) reshape(c2',3,1,[]) reshape(c2'.*0,3,1,[])];
infl = infl.*repmat(reshape(ori,1,1,[]),3,6,1);

end