 function [vort_edge] = fcnDVORTEDGE(idx, pts, theta, dvenum, valNELE, matROTANG, matEATT, matCENTER)

nedg = size(pts,1);

pts_one = fcnGLOBSTAR(pts(:,:,1) - matCENTER(dvenum(:,1),:), matROTANG(dvenum(:,1),:));
pts_two = fcnGLOBSTAR(pts(:,:,2) - matCENTER(dvenum(:,2),:), matROTANG(dvenum(:,2),:));

a2 = ones(nedg,2);
a1 = [pts_one(:,2) pts_two(:,2)];

b2 = ones(nedg,2);
b1 = [pts_one(:,1) pts_two(:,1)];

zer = zeros(nedg,1);

% A of element 1 carried into element 2 (A induces in xi, changes with eta)
dgamma1 = [a1(:,1), a2(:,1), zer, zer, zer];
dgamma2 = [a1(:,2), a2(:,2), zer, zer, zer].*-1;
vort_edge1 = fcnCREATEDSECT(sparse(nedg, valNELE*5), nedg, 5, matEATT(idx,1), matEATT(idx,2), dgamma1, dgamma2);

% A of element 1 carried into element 2 (A induces in xi, changes with eta)
dgamma1 = [zer, zer, b1(:,1), b2(:,1), zer];
dgamma2 = [zer, zer, b1(:,2), b2(:,2), zer].*-1;
vort_edge2 = fcnCREATEDSECT(sparse(nedg, valNELE*5), nedg, 5, matEATT(idx,1), matEATT(idx,2), dgamma1, dgamma2);

% % A of element 1 carried into element 2 (A induces in xi, changes with eta)
% dgamma1 = [zer, a2(:,1), zer, b2(:,1), zer];
% dgamma2 = [zer, a2(:,2), zer, b2(:,2), zer].*-1;
% vort_edge1 = fcnCREATEDSECT(sparse(nedg, valNELE*5), nedg, 5, dvenum(:,1), dvenum(:,2), dgamma1, dgamma2);
% 
% % A of element 1 carried into element 2 (A induces in xi, changes with eta)
% dgamma1 = [a1(:,1), zer, zer, zer, zer];
% dgamma2 = [a1(:,2).*cos(theta), zer, b1(:,2).*sin(theta), zer, zer].*-1;
% vort_edge2 = fcnCREATEDSECT(sparse(nedg, valNELE*5), nedg, 5, matEATT(idx,1), dvenum(:,2), dgamma1, dgamma2);
% 
% % B of element 1 carried into element 2 (B induces in eta, changes with xi)
% dgamma1 = [zer, zer, b1(:,1), zer, zer];
% dgamma2 = [a1(:,2).*-sin(theta), zer, b1(:,2).*cos(theta), zer, zer].*-1;
% vort_edge3 = fcnCREATEDSECT(sparse(nedg, valNELE*5), nedg, 5, matEATT(idx,1), dvenum(:,2), dgamma1, dgamma2);

vort_edge = [vort_edge1; vort_edge2];

end