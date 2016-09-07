function fpg = fcnTOGLOB(dvenum, fpl, DVE, DVECT, VLST)
% input (DVEnum, fp) in local and get back (fp) in global
% maybe works with velocity vectors?

fpg = fpl(:,1).*DVECT(dvenum,:,1) + fpl(:,2).*DVECT(dvenum,:,2) + fpl(:,3).*DVECT(dvenum,:,3) + VLST(DVE(dvenum,1,1),:);

end

