function [circ] = fcnDCIRC(idx, nedg, lambda, valNELE, matPLEX, matEATT, matELOC, matALIGN)
% Calculates circulation equations for the D-matrix

%(x,y) of all three vertices of HDVEs in local coordinates
x1 = reshape(matPLEX(1,1,matEATT(idx,:)),nedg,2);
x2 = reshape(matPLEX(2,1,matEATT(idx,:)),nedg,2);
x3 = reshape(matPLEX(3,1,matEATT(idx,:)),nedg,2);
y1 = reshape(matPLEX(1,2,matEATT(idx,:)),nedg,2);
y2 = reshape(matPLEX(2,2,matEATT(idx,:)),nedg,2);
y3 = reshape(matPLEX(3,2,matEATT(idx,:)),nedg,2);

lmb1 = reshape(lambda(matELOC(idx,:),1),nedg,2);
lmb2 = reshape(lambda(matELOC(idx,:),2),nedg,2);
lmb3 = reshape(lambda(matELOC(idx,:),3),nedg,2);

b2 = (lmb1.*x1+lmb2.*x2+lmb3.*x3);
b1 = b2.^2;
a2 = (lmb1.*y1+lmb2.*y2+lmb3.*y3);
a1 = a2.^2;
a3 = ones(nedg,2);
b3 = ones(nedg,2);

% align = abs(matALIGN(idx,:,:));
align = matALIGN(idx,:,:);

zer = zeros(length(a1(:,1)),1);

% gamma12 = [];
% gamma22 = [];
% gamma13 = [];
% gamma23 = [];
% 
% gamma12 = [a1(:,1), a2(:,1), a3(:,1), zer, zer, zer];
% gamma22 = [a1(:,2).*align(:,1,1), a2(:,2).*align(:,1,1), a3(:,2).*align(:,1,1), b1(:,2).*align(:,2,1), b2(:,2).*align(:,2,1), b3(:,2).*align(:,2,1)].*-1;
% 
% gamma13 = [zer, zer, zer, b1(:,1), b2(:,1), b3(:,1)];
% gamma23 = [a1(:,2).*align(:,1,2), a2(:,2).*align(:,1,2), a3(:,2).*align(:,1,2), b1(:,2).*align(:,2,2), b2(:,2).*align(:,2,2), b3(:,2).*align(:,2,2)].*-1;


% gamma22 = [a1(:,2).*align(:,1,1), a2(:,2).*align(:,1,1), a3(:,2).*align(:,1,1), b1(:,2).*align(:,2,1), b2(:,2).*align(:,2,1), b3(:,2).*align(:,2,1)].*-1;

gamma1 = [a1(:,1), a2(:,1), a3(:,1), b1(:,1), b2(:,1), b3(:,1)];
gamma2 = [a1(:,2), a2(:,2), a3(:,2), b1(:,2), b2(:,2), b3(:,2)].*-1;
circ = fcnCREATEDSECT(sparse(nedg, valNELE*6), nedg, 6, matEATT(idx,1), matEATT(idx,2), gamma1, gamma2);

  
end

