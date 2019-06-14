function [infl_glob] = fcnHDVEINDGLOB(dvenum, dvetype, fpg, matPLEX, matROTANG, matCONTROL, vecDVESYM, vecBOUNDIND)

infl_loc = fcnHDVEIND_DB(dvenum, dvetype, fpg, matPLEX, matROTANG, matCONTROL, vecBOUNDIND);

% Transforming
dvenum2 = reshape(repmat(dvenum,1,6,1)',[],1,1);
infl_tot = fcnSTARGLOB(reshape(permute(infl_loc,[2 3 1]),[],3,1), matROTANG(dvenum2,:));

% Symmetry
if any(vecDVESYM)
    idx = vecDVESYM(dvenum);
    infl_sym = fcnHDVEIND_DB(dvenum(idx), dvetype(idx), [fpg(idx,1) -fpg(idx,2) fpg(idx,3)], matPLEX, matROTANG, matCONTROL, vecBOUNDIND);

    % Transforming
    dvenum2 = reshape(repmat(dvenum(idx),1,6,1)',[],1,1);
    idx = vecDVESYM(dvenum2);
    infl_tot_sym = fcnSTARGLOB(reshape(permute(infl_sym,[2 3 1]),[],3,1), matROTANG(dvenum2,:));
    infl_tot_sym(:,2) = infl_tot_sym(:,2).*-1;
    infl_tot(idx,:) = infl_tot(idx,:) + infl_tot_sym;
end

infl_tot(isnan(infl_tot)) = 0;
infl_tot(isinf(infl_tot)) = 0;

infl_glob = reshape(infl_tot',3,6,[]);



end