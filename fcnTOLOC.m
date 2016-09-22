function fpl = fcnTOLOC(dvenum, fpg, DVE, DVECT, VLST)
% Given 3 global points and 3 local points, this transforms from global to local
% T.D.K 2016-09-05. 745-55 RIVER OAKS PLACE, SAN JOSE, CALIFORNIA, USA 95134

fpl(:,3) = dot(DVECT(dvenum,:,3), fpg - VLST(DVE(dvenum,1,1),:),2); % The projection of the vector to the fp onto the global normal of the element is the local zeta

r1 = (fpg - VLST(DVE(dvenum,1,1),:)) - (repmat(dot(fpg - VLST(DVE(dvenum,1,1),:), DVECT(dvenum,:,3),2),1,3)).*DVECT(dvenum,:,3); % Vector to FP projection onto plane of DVE from local origin

fpl(:,1) = dot(DVECT(dvenum,:,1), r1, 2); % eta
fpl(:,2) = dot(DVECT(dvenum,:,2), r1, 2); % xi

end