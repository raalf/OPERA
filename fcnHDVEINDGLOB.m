function [infl_glob] = fcnHDVEINDGLOB(dvenum, dvetype, fpg, matPLEX, matROTANG, matCONTROL)

infl_loc = fcnHDVEIND_DB(dvenum, dvetype, fpg, matPLEX, matROTANG, matCONTROL);

%% Transforming and Outputting
dvenum = reshape(repmat(dvenum,1,6,1)',[],1,1);
infl_tot = fcnSTARGLOB(reshape(permute(infl_loc,[2 3 1]),[],3,1), matROTANG(dvenum,:));

infl_tot(isnan(infl_tot)) = 0;
infl_tot(isinf(infl_tot)) = 0;

infl_glob = reshape(infl_tot',3,6,[]);

end