function [vort_edge] = fcnDVORTTE(idx, vnuma, nedg, lambda, valNELE, matPLEX, matEATT, e1vec)
% Calculates vorticity equations for the D-matrix

%(x,y) of all three vertices of HDVEs in local coordinates
x1 = reshape(matPLEX(1,1,nonzeros(matEATT(idx,:))),nedg,1);
x2 = reshape(matPLEX(2,1,nonzeros(matEATT(idx,:))),nedg,1);
x3 = reshape(matPLEX(3,1,nonzeros(matEATT(idx,:))),nedg,1);
y1 = reshape(matPLEX(1,2,nonzeros(matEATT(idx,:))),nedg,1);
y2 = reshape(matPLEX(2,2,nonzeros(matEATT(idx,:))),nedg,1);
y3 = reshape(matPLEX(3,2,nonzeros(matEATT(idx,:))),nedg,1);

% Barycentric coordinates (for vorticity, usually the vertices)
lmb1 = reshape(lambda([vnuma],1),nedg,1);
lmb2 = reshape(lambda([vnuma],2),nedg,1);
lmb3 = reshape(lambda([vnuma],3),nedg,1);

a2 = -ones(nedg,2);
a1 = -(lmb1.*y1+lmb2.*y2+lmb3.*y3);

b2 = ones(nedg,2);
b1 = (lmb1.*x1+lmb2.*x2+lmb3.*x3);

c2xi = (lmb1.*x1+lmb2.*x2+lmb3.*x3);
c2eta = (lmb1.*y1+lmb2.*y2+lmb3.*y3);

zer = a1(:,1).*0;

dgamma1 = [-a1(:,1).*e1vec(:,2), -a2(:,1).*e1vec(:,2), b1(:,1).*e1vec(:,1), b2(:,1).*e1vec(:,1), (-c2xi(:,1).*e1vec(:,2) + c2eta(:,1).*e1vec(:,1)), zer];

% Row indices of the rows where vorticity equations will go
rows = reshape([repmat([1:nedg]',1,6)]',[],1);

% Column indices for each circulation equation, col# = (DVE*6)-5 as each DVE gets a 6 column group
col1 = reshape([repmat([(nonzeros(matEATT(idx,:)).*6)-5],1,6)+repmat([0:5],nedg,1)]',[],1);

vort_edge = zeros(nedg, valNELE*6);
vort_edge(sub2ind(size(vort_edge),rows,col1)) = reshape(dgamma1',[],1);

end