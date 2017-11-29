function [curv] = fcnDCURV(idx, nedg, lambda, valNELE, matPLEX, matEATT, matELOC)
% Calculates curvature equations for the D-matrix


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

a2 = zeros(nedg,2);
a1 = ones(nedg,2);
b2 = zeros(nedg,2);
b1 = ones(nedg,2);
c3 = zeros(nedg,2);

ddgamma1 = [a1(:,1),a2(:,1),b1(:,1),b2(:,1),c3(:,1)];
ddgamma2 = [a1(:,2),a2(:,2),b1(:,2),b2(:,2),c3(:,2)].*-1;

% Row indices of the rows where circulation equations will go
rows = reshape([repmat([1:nedg]',1,5)]',[],1);

% Column indices for each circulation equation, col# = (DVE*5)-4 as each DVE gets a 6 column group
col1 = reshape([repmat([(matEATT(idx,1).*5)-4],1,5) + repmat([0:4],nedg,1)]',[],1);
col2 = reshape([repmat([(matEATT(idx,2).*5)-4],1,5) + repmat([0:4],nedg,1)]',[],1);

curv = zeros(nedg, valNELE*5);

curv(sub2ind(size(curv),rows,col1)) = reshape(ddgamma1',[],1);
curv(sub2ind(size(curv),rows,col2)) = reshape(ddgamma2',[],1);
    



end

