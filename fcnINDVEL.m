function [q_ind] = fcnINDVEL(dvenum, fpg, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, dvetype, matROTANG, matVSCOMB, matCENTER)
% This function takes in HDVE number and a corresponding global field point and returns an induced velocity
% in the global reference frame. 

% T.D.K 2016-09-11 6075 CUMULUS LANE, SAN DIEGO, CALIFORNIA, USA 92110

len = length(dvenum);

infl_glob = fcnHDVEINDFS(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX, dvetype, matROTANG, matVSCOMB, matCENTER);

q_ind = permute(sum(infl_glob.*repmat(reshape(matCOEFF(dvenum,:)',1,6,[]),3,1,1),2),[2 1 3]);

q_ind = reshape(permute(q_ind,[3 1 2]),[],3,1)./(4*pi);


end

