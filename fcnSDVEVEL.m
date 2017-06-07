function [q_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG, matVSCOMB, matCENTER)

len = length(fpg(:,1));
dvenum = reshape(repmat(1:valNELE,len,1),[],1);
dvetype = ones(size(dvenum));

fpg = repmat(fpg,valNELE,1);

[q_ind] = fcnINDVEL(dvenum, fpg, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, dvetype, matROTANG, matVSCOMB, matCENTER);

q_ind = reshape(sum(permute(reshape(q_ind',3,[],valNELE),[3 1 2]),1),3,[],1)';

end

