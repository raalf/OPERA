function [w_ind] = fcnWINDVEL(fpg, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER)
% This function finds the induced velocity at a point, due to the influence of the entire wake

len = length(fpg(:,1));
dvenum = reshape(repmat(1:valWNELE,len,1),[],1);
dvetype = ones(size(dvenum)) + 1;

oldest_wake = [1:valWSIZE] + valWSIZE;
dvetype(oldest_wake) = 3;

fpg = repmat(fpg,valWNELE,1);

[q_ind] = fcnINDVEL(dvenum, dvetype, fpg, matWCOEFF, matWPLEX, matWROTANG, matWCENTER);

w_ind = reshape(sum(reshape(q_ind', len*3, [])',1),3,[])';

end

