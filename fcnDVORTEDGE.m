 function [vort_edge] = fcnDVORTEDGE(pts, dvenum, valNELE, matROTANG, matCENTER)

nedg = size(pts,1);

pts_one = fcnGLOBSTAR(pts(:,:,1) - matCENTER(dvenum(:,1),:), matROTANG(dvenum(:,1),:));
pts_two = fcnGLOBSTAR(pts(:,:,2) - matCENTER(dvenum(:,2),:), matROTANG(dvenum(:,2),:));

a2 = ones(nedg,2);
a1 = [pts_one(:,2) pts_two(:,2)];
c2a = [pts_one(:,1) pts_two(:,1)];

b2 = ones(nedg,2);
b1 = [pts_one(:,1) pts_two(:,1)];
c2b = [pts_one(:,2) pts_two(:,2)];

zer = zeros(nedg,1);

% A of element 1 carried into element 2 (A induces in xi, changes with eta)
dgamma1 = [a1(:,1), a2(:,1), zer, zer, c2a(:,1), zer];
dgamma2 = [a1(:,2), a2(:,2), zer, zer, c2a(:,2), zer].*-1;
vort_edge1 = fcnCREATEDSECT(sparse(nedg, valNELE*6), nedg, 6, dvenum(:,1), dvenum(:,2), dgamma1, dgamma2);

% A of element 1 carried into element 2 (A induces in xi, changes with eta)
dgamma1 = [zer, zer, b1(:,1), b2(:,1), c2b(:,1), zer];
dgamma2 = [zer, zer, b1(:,2), b2(:,2), c2b(:,2), zer].*-1;
vort_edge2 = fcnCREATEDSECT(sparse(nedg, valNELE*6), nedg, 6, dvenum(:,1), dvenum(:,2), dgamma1, dgamma2);

vort_edge = [vort_edge1; vort_edge2];

end