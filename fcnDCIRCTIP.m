function [circ_tip] = fcnDCIRCTIP(idx, nedg, lambda_vert, lambda_mid, valNELE, matPLEX, matEATT, matELOC, vnumc)
% Calculates circulation equations for the D-matrix

% (x,y) of all three vertices of HDVEs in local coordinates
x1 = repmat(reshape(matPLEX(1,1,nonzeros(matEATT(idx,:))),nedg,1),3,1);
x2 = repmat(reshape(matPLEX(2,1,nonzeros(matEATT(idx,:))),nedg,1),3,1);
x3 = repmat(reshape(matPLEX(3,1,nonzeros(matEATT(idx,:))),nedg,1),3,1);
y1 = repmat(reshape(matPLEX(1,2,nonzeros(matEATT(idx,:))),nedg,1),3,1);
y2 = repmat(reshape(matPLEX(2,2,nonzeros(matEATT(idx,:))),nedg,1),3,1);
y3 = repmat(reshape(matPLEX(3,2,nonzeros(matEATT(idx,:))),nedg,1),3,1);

% Barycentric coordinates (the edge midpoints)
lmb1 = reshape(lambda_mid(nonzeros(matELOC(idx,:)),1),nedg,1);
lmb2 = reshape(lambda_mid(nonzeros(matELOC(idx,:)),2),nedg,1);
lmb3 = reshape(lambda_mid(nonzeros(matELOC(idx,:)),3),nedg,1);

% Barycentric coordinates (the vertices)
lmb1 = [lmb1; lambda_vert(vnumc(:,1),1); lambda_vert(vnumc(:,2),1)];
lmb2 = [lmb2; lambda_vert(vnumc(:,1),2); lambda_vert(vnumc(:,2),2)];
lmb3 = [lmb3; lambda_vert(vnumc(:,1),3); lambda_vert(vnumc(:,2),3)];

a2 = (lmb1.*x1+lmb2.*x2+lmb3.*x3);
a1 = a2.^2;
b2 = (lmb1.*y1+lmb2.*y2+lmb3.*y3);
b1 = b2.^2;
a3 = ones(nedg*3,1);
b3 = ones(nedg*3,1);

gamma_tip = [a1(:,1), a2(:,1), a3(:,1), b1(:,1), b2(:,1), b3(:,1)];

rows = reshape([repmat([1:nedg*3]',1,6)]',[],1);
col4 = reshape([repmat([(nonzeros(matEATT(idx,:)).*6)-5],3,6)+repmat([0:5],nedg*3,1)]',[],1);

circ_tip = zeros(nedg*3, valNELE*6);
circ_tip(sub2ind(size(circ_tip),rows,col4)) = reshape(gamma_tip',[],1);


end

