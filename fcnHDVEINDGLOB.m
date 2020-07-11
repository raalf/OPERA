function [infl_glob] = fcnHDVEINDGLOB(dvenum, dvetype, fpg, matPLEX, matROTANG, matCONTROL, vecBOUNDIND, ztol)

infl_loc = fcnHDVEIND(dvenum, dvetype, fpg, matPLEX, matROTANG, matCONTROL, vecBOUNDIND, ztol);

% Transforming
dvenum2 = reshape(repmat(dvenum,1,6,1)',[],1,1);
infl_tot = fcnSTARGLOB(reshape(permute(infl_loc,[2 3 1]),[],3,1), matROTANG(dvenum2,:));

infl_tot(isnan(infl_tot)) = 0;
infl_tot(isinf(infl_tot)) = 0;

infl_glob = reshape(infl_tot',3,6,[]);

end