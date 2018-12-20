function [q_ind] = fcnINDVEL(dvenum, dvetype, fpg, matCOEFF, matPLEX, matROTANG, matCENTER)
% This function takes in HDVE number and a corresponding global field point and returns an induced velocity
% in the global reference frame. 
% T.D.K 2016-09-11 6075 CUMULUS LANE, SAN DIEGO, CALIFORNIA, USA 92110

% infl_loc = fcnHDVEIND(dvenum, dvetype, fpg, matPLEX, matROTANG, matCENTER);
% infl_loc = fcnHDVEIND_VS(dvenum, dvetype, fpg, matPLEX, matROTANG, matCENTER);
infl_loc = fcnHDVEIND_DB(dvenum, dvetype, fpg, matPLEX, matROTANG, matCENTER);

q_ind = permute(sum(infl_loc.*repmat(reshape(matCOEFF(dvenum,:)',1,5,[]),3,1,1),2),[2 1 3]);
q_ind = reshape(permute(q_ind,[3 1 2]),[],3,1)./(-4*pi);

q_ind = fcnSTARGLOB(q_ind, matROTANG(dvenum,:));
end

