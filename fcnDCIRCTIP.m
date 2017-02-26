function [circ_tip] = fcnDCIRCTIP(idx, nedg, lambda, valNELE, matPLEX, matEATT, matELOC)
% Calculates circulation equations for the D-matrix

% (x,y) of all three vertices of HDVEs in local coordinates
x1 = reshape(matPLEX(1,1,nonzeros(matEATT(idx,:))),nedg,1);
x2 = reshape(matPLEX(2,1,nonzeros(matEATT(idx,:))),nedg,1);
x3 = reshape(matPLEX(3,1,nonzeros(matEATT(idx,:))),nedg,1);
y1 = reshape(matPLEX(1,2,nonzeros(matEATT(idx,:))),nedg,1);
y2 = reshape(matPLEX(2,2,nonzeros(matEATT(idx,:))),nedg,1);
y3 = reshape(matPLEX(3,2,nonzeros(matEATT(idx,:))),nedg,1);

lmb1 = reshape(lambda(nonzeros(matELOC(idx,:)),1),nedg,1);
lmb2 = reshape(lambda(nonzeros(matELOC(idx,:)),2),nedg,1);
lmb3 = reshape(lambda(nonzeros(matELOC(idx,:)),3),nedg,1);

a2 = (lmb1.*x1+lmb2.*x2+lmb3.*x3);
a1 = a2.^2;
b2 = (lmb1.*y1+lmb2.*y2+lmb3.*y3);
b1 = b2.^2;

c3 = ones(nedg,2);

gamma_tip = [a1(:,1), a2(:,1), b1(:,1), b2(:,1), c3(:,1)];

rows = reshape([repmat([1:nedg]',1,5)]',[],1);
col4 = reshape([repmat([(nonzeros(matEATT(idx,:)).*5)-4],1,5)+repmat([0:4],nedg,1)]',[],1);

circ_tip = zeros(nedg, valNELE*5);
circ_tip(sub2ind(size(circ_tip),rows,col4)) = reshape(gamma_tip',[],1);


end

