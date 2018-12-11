function [vort_edge] = fcnDVORTTE(idx, pts, theta, dvenum, valNELE, matROTANG, matEATT, matCENTER)

nedg = size(pts,1);

pts_one = fcnGLOBSTAR(pts(:,:,1) - matCENTER(dvenum(:,1),:), matROTANG(dvenum(:,1),:));

a2 = ones(nedg,2);
a1 = pts_one(:,2);

b2 = ones(nedg,2);
b1 = pts_one(:,1);

zer = zeros(nedg,1);

dgamma1 = [zer, a2(:,1).*cos(theta), zer, b2(:,1).*sin(theta), zer];
vort_edge1 = fcnCREATEDSECT(sparse(nedg, valNELE*5), nedg, 5, dvenum, [], dgamma1, []);

dgamma1 = [a1(:,1).*cos(theta), zer, b1(:,1).*sin(theta), zer, zer];
vort_edge2 = fcnCREATEDSECT(sparse(nedg, valNELE*5), nedg, 5, dvenum, [], dgamma1, []);

vort_edge = real([vort_edge1; vort_edge2]);

end