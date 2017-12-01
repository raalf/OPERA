function [vort] = fcnDVORTVERT(idx, vnuma, vnumb, nedg, lambda, valNELE, matPLEX, matEATT, matELOC, matALIGN)
% Calculates vorticity equations for the D-matrix

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

a2 = ones(nedg,2);
a1 = 2.*(lmb1.*y1+lmb2.*y2+lmb3.*y3);

b2 = ones(nedg,2);
b1 = 2.*(lmb1.*x1+lmb2.*x2+lmb3.*x3);

c2eta = (lmb1.*y1+lmb2.*y2+lmb3.*y3);
c2xi = (lmb1.*x1+lmb2.*x2+lmb3.*x3);

% align = abs(matALIGN(idx,:,:));
align = matALIGN(idx,:,:);

zer = a1(:,1).*0;
dgamma12 = [a1(:,1), a2(:,1), zer, zer, c2xi(:,1), zer];
dgamma22 = [a1(:,2), a2(:,2), zer, zer, c2xi(:,2), zer].*-1;

dgamma13 = [zer, zer, b1(:,1), b2(:,1), c2eta(:,1), zer];
dgamma23 = [zer, zer, b1(:,2), b2(:,2), c2eta(:,2), zer].*-1;


lines = nedg*2;

% % Row indices of the rows where circulation equations will go
rows = reshape([repmat([1:lines]',1,6)]',[],1);

% Column indices for each circulation equation, col# = (DVE*5)-4 as each DVE gets a 6 column group
col1 = reshape([repmat([(matEATT(idx,1).*6)-5],2,6) + repmat([0:5],lines,1)]',[],1);
col2 = reshape([repmat([(matEATT(idx,2).*6)-5],2,6) + repmat([0:5],lines,1)]',[],1);

vort = zeros(lines, valNELE*6);

dgamma1 = [dgamma12; dgamma13];
dgamma2 = [dgamma22; dgamma23];
vort(sub2ind(size(vort),rows,col1)) = reshape(dgamma1',[],1);
vort(sub2ind(size(vort),rows,col2)) = reshape(dgamma2',[],1);

end