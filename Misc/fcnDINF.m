function [dinf] = fcnDINF(idx, nedg, dve1, dve2, fpg, matDVE, matDVECT, matVLST, matPLEX, dvetype, matROTANG, matEATT, valNELE)

len = length(fpg(:,1));

[a1, a2, a3, b1, b2, b3] = fcnHDVEIND(dve1, fpg, matDVE, matDVECT, matVLST, matPLEX, dvetype, matROTANG);

temp1 = [a1 a2 a3 b1 b2 b3];
temp1 = reshape(reshape(temp1', 1, 18, []), 3, 6, len);
temp1 = reshape(permute(temp1, [2 1 3]), size(temp1, 2), [])';

[a1, a2, a3, b1, b2, b3] = fcnHDVEIND(dve2, fpg, matDVE, matDVECT, matVLST, matPLEX, dvetype, matROTANG);

temp2 = [a1 a2 a3 b1 b2 b3].*-1;
temp2 = reshape(reshape(temp2', 1, 18, []), 3, 6, len);
temp2 = reshape(permute(temp1, [2 1 3]), size(temp2, 2), [])';

% % Row indices of the rows where circulation equations will go
rows = reshape([repmat([1:nedg*3]',1,6)]',[],1);

% Column indices for each circulation equation, col# = (DVE*5)-4 as each DVE gets a 6 column group
col1 = reshape([repmat([(matEATT(idx,1).*6)-5],1,18) + repmat([0:5],nedg*1,3)]',[],1);
col2 = reshape([repmat([(matEATT(idx,2).*6)-5],1,18) + repmat([0:5],nedg*1,3)]',[],1);

dinf = zeros(nedg*3, valNELE*6);

dinf(sub2ind(size(dinf),rows,col1)) = reshape(temp1',[],1);
dinf(sub2ind(size(dinf),rows,col2)) = reshape(temp2',[],1);

end

