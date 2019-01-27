function [circ_tip] = fcnDCIRC2(pts, dvenum, valNELE, matROTANG, matCENTER)

nedg = size(pts,1);

pts_one = fcnGLOBSTAR(pts(:,:,1) - matCENTER(dvenum(:,1),:), matROTANG(dvenum(:,1),:));

a2 = [pts_one(:,2)];
a1 = 0.5.*(a2.^2);

b2 = [pts_one(:,1)];
b1 = 0.5.*(b2.^2);

c2 = [pts_one(:,1).*pts_one(:,2)];
c3 = ones(nedg,2);

gamma1 = [a1(:,1), a2(:,1), b1(:,1), b2(:,1), c2(:,1), c3(:,1)];

circ_tip = fcnCREATEDSECT(sparse(nedg, valNELE*6), nedg, 6, dvenum, [], gamma1, []);
end

