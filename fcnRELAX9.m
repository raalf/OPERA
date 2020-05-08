function [matWELST, matWVLST, matWPLEX, matWDVECT, matWCENTER, matWROTANG] = ...
    fcnRELAX9(valTIMESTEP, valDELTIME, valNELE, matCOEFF, matPLEX, valWNELE, matWCOEFF, matWDVE, matWVLST, matWPLEX, valWSIZE, ...
    matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST, matWEIDX, vecDVESYM, vecWDVESYM, matWVGRID, vecWDVEFLIP, ...
    valPRESTEPS, matWDVEGRID, valWOFF, vecWOFF, vecWNORM, vecHUB, matVLST, matDVE, vecDVEFLIP)

%% Finding induced velocities at all wake vertices
sz = 10;
q_ind = matWCENTER.*0;
for i = 1:sz:valTIMESTEP
    if i + sz*2 > valTIMESTEP
        idx_rows = [i:valTIMESTEP]';
    else
        idx_rows = [i:i+sz]';
    end
    
    q_ind(matWDVEGRID(idx_rows,:),:) = fcnSDVEVEL(matWCENTER(matWDVEGRID(idx_rows,:),:), valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM, [], 0) + fcnSDVEVEL(matWCENTER(matWDVEGRID(idx_rows,:),:), valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, vecWDVESYM, [], 0);
    
    if ~isnan(valWOFF)
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
        
%         hold on
%         scatter3(tmpVLST(:,1), tmpVLST(:,2), tmpVLST(:,3), 'sk');
%         scatter3(tmpWVLST(:,1), tmpWVLST(:,2), tmpWVLST(:,3), 'ob');
%         quiver3(tmpCENTER(:,1), tmpCENTER(:,2), tmpCENTER(:,3), DNORM(:,1), DNORM(:,2), DNORM(:,3),'r');
%         quiver3(tmpWCENTER(:,1), tmpWCENTER(:,2), tmpWCENTER(:,3), WDNORM(:,1), WDNORM(:,2), WDNORM(:,3),'b');
%         hold off
        
        
        q_ind(matWDVEGRID(idx_rows,:),:) = q_ind(matWDVEGRID(idx_rows,:),:) + ...
            fcnSDVEVEL(matWCENTER(matWDVEGRID(idx_rows,:),:), valNELE, matCOEFF, matPLEX, tmpROTANG, tmpCENTER, [], [], 0) + ...
            fcnSDVEVEL(matWCENTER(matWDVEGRID(idx_rows,:),:), valWNELE, matWCOEFF, matWPLEX, tmpWROTANG, tmpWCENTER, [], [], 0);
        
    end
    
    
end

%% Moving wake vertices
% matWCENTER = matWCENTER + q_ind.*valDELTIME;
tmp2 = matWVLST.*0;
vcount = zeros(size(matWVLST,1),1);
total_weight = matWVLST(:,1).*0;
for i = 1:3
    dist = sqrt(sum((matWVLST(matWDVE(:,i),:) - matWCENTER).^2,2));
    tmp2(matWDVE(:,i),:) = tmp2(matWDVE(:,i),:) + q_ind.*dist;
    [gc,gr] = groupcounts(matWDVE(:,i));
    vcount(gr,:) = vcount(gr,:) + gc;
    total_weight(matWDVE(:,i),:) = total_weight(matWDVE(:,i),:) + dist;
end
tmp2 = tmp2./(vcount.*total_weight);
tmp2(isnan(tmp2)) = 0;

tmp2(matWELST(vecWLE,:),:) = tmp2(matWELST(vecWLE,:),:).*0; % Trailing edge of wing
tmp2(matWELST(vecWTE,:),:) = tmp2(matWELST(vecWTE,:),:).*0; % Trailing edge of LE wake element row
tmp2(matWVGRID(end,:),:) = tmp2(matWVGRID(end-1,:),:); % Velocities of oldest row of wake vertices are set equal to second oldest
matWVLST = matWVLST + tmp2.*valDELTIME;

%% Recreating wake point matrix, and regenerating wake HDVE parameters
% matWCENTER
matWCENTER = (matWVLST(matWDVE(:,1),:) + matWVLST(matWDVE(:,2),:) + matWVLST(matWDVE(:,3),:))./3;

% matWPLEX, matWDVECT, matWROTANG
P = permute(reshape(matWVLST(matWDVE(:,:)',:)', 3, 3, []), [2 1 3]);
DNORM = cross(matWVLST(matWDVE(:,2),:) - matWVLST(matWDVE(:,3),:), matWVLST(matWDVE(:,1),:) - matWVLST(matWDVE(:,3),:), 2);
DNORM = DNORM./sqrt(sum(DNORM.^2,2));
DNORM(vecWDVEFLIP,:) = DNORM(vecWDVEFLIP,:).*-1;
[matWPLEX, matWDVECT, matWROTANG] = fcnTRITOLEX(P, DNORM, matWCENTER);

end
