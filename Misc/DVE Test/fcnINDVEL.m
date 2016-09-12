function [q_ind] = fcnINDVEL(dvenum, fpg, COEFF, DVE, DVECT, VLST, DNORM, PLEX)
% This function takes in HDVE number and a corresponding global field point and returns an induced velocity
% in the global reference frame. 

% T.D.K 2016-09-11 6075 CUMULUS LANE, SAN DIEGO, CALIFORNIA, USA 92110

len = length(dvenum);

[a1, a2, b1, b2, c3] = fcnHDVEIND(dvenum, fpg, DVE, DVECT, VLST, DNORM, PLEX);

D = [a1 a2 b1 b2 c3];
D = reshape(reshape(D', 1, 15, []), 3, 5, len);

q_ind = permute(sum(D.*COEFF,2),[2 1 3]);
q_ind = reshape(permute(q_ind,[3 1 2]),[],3,1);

end

