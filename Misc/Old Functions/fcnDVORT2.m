function [vort_e, vort_x] = fcnDVORT2(idx, nedg, lambda, valNELE, matPLEX, matEATT, matELOC, matALIGN)
% Calculates vorticity equations for the D-matrix

%(x,y) of all three vertices of HDVEs in local coordinates
x1 = reshape(matPLEX(1,1,matEATT(idx,:)),nedg,2);
x2 = reshape(matPLEX(2,1,matEATT(idx,:)),nedg,2);
x3 = reshape(matPLEX(3,1,matEATT(idx,:)),nedg,2);
y1 = reshape(matPLEX(1,2,matEATT(idx,:)),nedg,2);
y2 = reshape(matPLEX(2,2,matEATT(idx,:)),nedg,2);
y3 = reshape(matPLEX(3,2,matEATT(idx,:)),nedg,2);

% Barycentric coordinates (for vorticity here, edge midpoints)
lmb1 = reshape(lambda(matELOC(idx,:),1),nedg,2);
lmb2 = reshape(lambda(matELOC(idx,:),2),nedg,2);
lmb3 = reshape(lambda(matELOC(idx,:),3),nedg,2);

a2 = ones(nedg,2);
a1 = 2.*(lmb1.*x1+lmb2.*x2+lmb3.*x3);
b2 = ones(nedg,2);
b1 = 2.*(lmb1.*y1+lmb2.*y2+lmb3.*y3);
c3 = zeros(nedg,2);

% Eta direction
dgamma1 = [a1(:,1), a2(:,1), zeros(nedg,1), zeros(nedg,1), zeros(nedg,1)].*-1;
dgamma2 = [a1(:,2).*matALIGN(idx,1,1), a2(:,2).*matALIGN(idx,1,1), b1(:,2).*matALIGN(idx,2,1), b2(:,2).*matALIGN(idx,2,1), zeros(nedg,1)];

% Row indices of the rows where vorticity equations will go
rows = reshape([repmat([1:nedg]',1,5)]',[],1);

% Column indices for each circulation equation, col# = (DVE*6)-5 as each DVE gets a 6 column group
col1 = reshape([repmat([(matEATT(idx,1).*5)-4],1,5)+repmat([0:4],nedg,1)]',[],1);
col2 = reshape([repmat([(matEATT(idx,2).*5)-4],1,5)+repmat([0:4],nedg,1)]',[],1);

vort_e = zeros(nedg, valNELE*5);
vort_e(sub2ind(size(vort_e),rows,col1)) = reshape(dgamma1',[],1);
vort_e(sub2ind(size(vort_e),rows,col2)) = reshape(dgamma2',[],1);

% Xsi Direction
dgamma1 = [zeros(nedg,1), zeros(nedg,1), b1(:,1), b2(:,1), zeros(nedg,1)].*-1;
dgamma2 = [a1(:,2).*matALIGN(idx,1,2), a2(:,2).*matALIGN(idx,1,2), b1(:,2).*matALIGN(idx,2,2), b2(:,2).*matALIGN(idx,2,2), zeros(nedg,1)];

% Row indices of the rows where vorticity equations will go
rows = reshape([repmat([1:nedg]',1,5)]',[],1);

% Column indices for each circulation equation, col# = (DVE*6)-5 as each DVE gets a 6 column group
col1 = reshape([repmat([(matEATT(idx,1).*5)-4],1,5)+repmat([0:4],nedg,1)]',[],1);
col2 = reshape([repmat([(matEATT(idx,2).*5)-4],1,5)+repmat([0:4],nedg,1)]',[],1);

vort_x = zeros(nedg, valNELE*5);
vort_x(sub2ind(size(vort_x),rows,col1)) = reshape(dgamma1',[],1);
vort_x(sub2ind(size(vort_x),rows,col2)) = reshape(dgamma2',[],1);

end