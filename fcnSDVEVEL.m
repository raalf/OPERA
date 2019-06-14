function [q_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM, vecBOUNDIND)
% T.D.K 2019-01-12 AC1886 OVER IOWA @ 35,000 FT

len = length(fpg(:,1));
dvenum = reshape(repmat(1:valNELE,len,1),[],1);
dvetype = ones(size(dvenum));

fpg = repmat(fpg,valNELE,1);

if ~isempty(vecBOUNDIND)
    tmp = vecBOUNDIND(dvenum);
else
    tmp = [];
end
[q_ind] = fcnINDVEL(dvenum, dvetype, fpg, matCOEFF, matPLEX, matROTANG, matCENTER, [], tmp);

if any(vecDVESYM)
    idx = vecDVESYM(dvenum);
    q_sym = fcnINDVEL(dvenum(idx), dvetype(idx), [fpg(idx,1) -fpg(idx,2) fpg(idx,3)], matCOEFF, matPLEX, matROTANG, matCENTER, [], vecBOUNDIND(idx));
    q_sym = q_sym.*[1 -1 1];
    q_ind(idx,:) = q_ind(idx,:) + q_sym;
end

q_ind = reshape(sum(permute(reshape(q_ind',3,[],valNELE),[3 1 2]),1),3,[],1)';

end

