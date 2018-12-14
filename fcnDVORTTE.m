function [vort_edge] = fcnDVORTTE(idx, pts, dvenum, valNELE, matROTANG, matEATT, matCENTER)

nedg = size(pts,1);

% pts_one = fcnGLOBSTAR(pts(:,:,1) - matCENTER(dvenum(:,1),:), matROTANG(dvenum(:,1),:));

a2 = ones(nedg,2);
% a1 = pts_one(:,2);
a1 = a2;
b2 = ones(nedg,2);
% b1 = pts_one(:,1);
b1 = b2;
zer = zeros(nedg,1);

dgamma1 = [zer, zer zer, b2(:,1), zer];
vort_edge1 = fcnCREATEDSECT(sparse(nedg, valNELE*5), nedg, 5, dvenum, [], dgamma1, []);

dgamma1 = [zer, zer, b1(:,1), zer, zer];
vort_edge2 = fcnCREATEDSECT(sparse(nedg, valNELE*5), nedg, 5, dvenum, [], dgamma1, []);

% dgamma1 = [zer, a2(:,1), zer, zer, zer];
% vort_edge1 = fcnCREATEDSECT(sparse(nedg, valNELE*5), nedg, 5, dvenum, [], dgamma1, []);
% 
% dgamma1 = [a1(:,1), zer, zer, zer, zer];
% vort_edge2 = fcnCREATEDSECT(sparse(nedg, valNELE*5), nedg, 5, dvenum, [], dgamma1, []);

vort_edge = real([vort_edge1; vort_edge2]);

end