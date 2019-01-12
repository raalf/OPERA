function [q_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER)
% T.D.K 2019-01-12 AC1886 OVER IOWA @ 35,000 FT

len = length(fpg(:,1));
dvenum = reshape(repmat(1:valNELE,len,1),[],1);
dvetype = ones(size(dvenum));

fpg = repmat(fpg,valNELE,1);

[q_ind] = fcnINDVEL(dvenum, dvetype, fpg, matCOEFF, matPLEX, matROTANG, matCENTER);
q_ind = reshape(sum(permute(reshape(q_ind',3,[],valNELE),[3 1 2]),1),3,[],1)';

end

