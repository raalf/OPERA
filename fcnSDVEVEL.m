function [q_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX)

len = length(fpg(:,1));
dvenum = repmat([1:valNELE]',len,1);

fpg = repmat(fpg,valNELE,1);

[q_ind] = fcnINDVEL(dvenum, fpg, matCOEFF, matDVE, matDVECT, matVLST, matPLEX);

q_ind = reshape(sum(permute(reshape(q_ind',3,[],valNELE),[2 1 3]),1),3,[],1)';

end

