function fpg = fcnTOGLOB(dvenum, fpl, DVE, DVECT, VLST)
% input (DVEnum, fp) in local and get back (fp) in global
% maybe works with velocity vectors?
% T.D.K 2016-09-05. 745-55 RIVER OAKS PLACE, SAN JOSE, CALIFORNIA, USA 95134

fpg = fpl(:,1).*DVECT(dvenum,:,1) + fpl(:,2).*DVECT(dvenum,:,2) + fpl(:,3).*DVECT(dvenum,:,3) + VLST(DVE(dvenum,1,1),:);

end

