function fpl = fcnTOLOC(dvenum, fpg, DVE, DVECT, VLST, DNORM)
% Given 3 global points and 3 local points, I need to create a transformation from global to local

fpl(:,3) = dot(DNORM(dvenum,:), fpg - VLST(DVE(dvenum,1,1),:),2); % The projection of the vector to the fp onto the global normal of the element is the local zeta

r1 = (fpg - fpl(:,3).*DNORM(dvenum,:)) - VLST(DVE(dvenum,1,1),:); % Vector to FP projection onto plane of DVE from local origin

fpl(:,1) = dot(DVECT(dvenum,:,1), r1, 2); % eta
fpl(:,2) = dot(DVECT(dvenum,:,2), r1, 2); % xi

end