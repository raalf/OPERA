function [tmpROTANG, tmpCENTER, tmpWROTANG, tmpWCENTER] = ...
    fcnMIRROR(matVLST, matDVE, vecDVEFLIP, matWVLST, matWDVE, vecWDVEFLIP, ...
    vecHUB, valWOFF, vecWOFF, vecWNORM)

tmpVLST = fcnMIRRORPTS(matVLST, vecHUB, valWOFF, vecWOFF, vecWNORM);
tmpCENTER = (tmpVLST(matDVE(:,1),:) + tmpVLST(matDVE(:,2),:) + tmpVLST(matDVE(:,3),:))./3;
P = permute(reshape(tmpVLST(matDVE(:,:)',:)', 3, 3, []), [2 1 3]);
DNORM = cross(tmpVLST(matDVE(:,2),:) - tmpVLST(matDVE(:,3),:), tmpVLST(matDVE(:,1),:) - tmpVLST(matDVE(:,3),:), 2);
DNORM = DNORM./sqrt(sum(DNORM.^2,2));
DNORM(~vecDVEFLIP,:) = DNORM(~vecDVEFLIP,:).*-1;
[~, ~, tmpROTANG] = fcnTRITOLEX(P, DNORM, tmpCENTER);

tmpWVLST = fcnMIRRORPTS(matWVLST, vecHUB, valWOFF, vecWOFF, vecWNORM);
tmpWCENTER = (tmpWVLST(matWDVE(:,1),:) + tmpWVLST(matWDVE(:,2),:) + tmpWVLST(matWDVE(:,3),:))./3;
P = permute(reshape(tmpWVLST(matWDVE(:,:)',:)', 3, 3, []), [2 1 3]);
WDNORM = cross(tmpWVLST(matWDVE(:,2),:) - tmpWVLST(matWDVE(:,3),:), tmpWVLST(matWDVE(:,1),:) - tmpWVLST(matWDVE(:,3),:), 2);
WDNORM = WDNORM./sqrt(sum(WDNORM.^2,2));
WDNORM(~vecWDVEFLIP,:) = WDNORM(~vecWDVEFLIP,:).*-1;
[~, ~, tmpWROTANG] = fcnTRITOLEX(P, WDNORM, tmpWCENTER);

% hold on
% scatter3(tmpVLST(:,1), tmpVLST(:,2), tmpVLST(:,3), 'sk');
% scatter3(tmpWVLST(:,1), tmpWVLST(:,2), tmpWVLST(:,3), 'ob');
% quiver3(tmpCENTER(:,1), tmpCENTER(:,2), tmpCENTER(:,3), DNORM(:,1), DNORM(:,2), DNORM(:,3),'r');
% quiver3(tmpWCENTER(:,1), tmpWCENTER(:,2), tmpWCENTER(:,3), WDNORM(:,1), WDNORM(:,2), WDNORM(:,3),'b');
% hold off

end