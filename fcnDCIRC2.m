function [circ] = fcnDCIRC2(idx, vnuma, vnumb, nedg, lambda, valNELE, matPLEX, matEATT, matELOC)
% Calculates circulation equations for the D-matrix

% A1 A2 B1 B2 C2 C3
% A changes with eta, B changes with xi, C changes with both

%(x,y) of all three vertices of HDVEs in local coordinates
x1 = reshape(matPLEX(1,1,matEATT(idx,:)),nedg,2);
x2 = reshape(matPLEX(2,1,matEATT(idx,:)),nedg,2);
x3 = reshape(matPLEX(3,1,matEATT(idx,:)),nedg,2);
y1 = reshape(matPLEX(1,2,matEATT(idx,:)),nedg,2);
y2 = reshape(matPLEX(2,2,matEATT(idx,:)),nedg,2);
y3 = reshape(matPLEX(3,2,matEATT(idx,:)),nedg,2);

% Barycentric coordinates (for vorticity, usually the vertices)
lmb1 = reshape(lambda([vnuma vnumb],1),nedg,2);
lmb2 = reshape(lambda([vnuma vnumb],2),nedg,2);
lmb3 = reshape(lambda([vnuma vnumb],3),nedg,2);

a2 = (lmb1.*y1+lmb2.*y2+lmb3.*y3);
a1 = 0.5.*(a2.^2);

b2 = (lmb1.*x1+lmb2.*x2+lmb3.*x3);
b1 = 0.5.*(b2.^2);

c2 = a2.*b2;
c3 = ones(nedg,2);

gamma1 = [a1(:,1), a2(:,1), b1(:,1), b2(:,1), c2(:,1), c3(:,1)];
gamma2 = [a1(:,2), a2(:,2), b1(:,2), b2(:,2), c2(:,2), c3(:,2)].*-1;

circ = fcnCREATEDSECT(sparse(nedg, valNELE*6), nedg, 6, matEATT(idx,1), matEATT(idx,2), gamma1, gamma2);

end

