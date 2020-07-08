function [matWELST, matWVLST, matWPLEX, matWDVECT, matWCENTER, matWROTANG] = ...
    fcnRELAX9(valTIMESTEP, valDELTIME, valNELE, matCOEFF, matPLEX, valWNELE, matWCOEFF, ... 
    matWDVE, matWVLST, matWPLEX, matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST, ...
    vecDVESYM, vecWDVESYM, matWVGRID, vecWDVEFLIP, matWDVEGRID)

%% Finding induced velocities at all wake vertices
sz = 10;
q_ind = matWCENTER.*0;
for i = 1:sz:valTIMESTEP
    if i + sz*2 > valTIMESTEP
        idx_rows = [i:valTIMESTEP]';
    else
        idx_rows = [i:i+sz]';
    end
    
%     q_ind(matWDVEGRID(idx_rows,:),:) = (fcnSDVEVEL(matWCENTER(matWDVEGRID(idx_rows,:),:) + matWDVECT(matWDVEGRID(idx_rows,:),:,3).*valZTOL, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM, [], 0) + ...
%         fcnSDVEVEL(matWCENTER(matWDVEGRID(idx_rows,:),:) - matWDVECT(matWDVEGRID(idx_rows,:),:,3).*valZTOL, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM, [], 0))./2 + ...        
%         (fcnSDVEVEL(matWCENTER(matWDVEGRID(idx_rows,:),:) + matWDVECT(matWDVEGRID(idx_rows,:),:,3).*valZTOL, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, vecWDVESYM, [], 0) + ...
%         fcnSDVEVEL(matWCENTER(matWDVEGRID(idx_rows,:),:) - matWDVECT(matWDVEGRID(idx_rows,:),:,3).*valZTOL, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, vecWDVESYM, [], 0))./2;
 
    q_ind(matWDVEGRID(idx_rows,:),:) = fcnSDVEVEL(matWCENTER(matWDVEGRID(idx_rows,:),:), valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM, [], 0) + ...
        fcnSDVEVEL(matWCENTER(matWDVEGRID(idx_rows,:),:), valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, vecWDVESYM, [], 0);    
end

hold on
quiver3(matWCENTER(:,1), matWCENTER(:,2), matWCENTER(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3));
hold off

%% Moving wake vertices
% matWCENTER = matWCENTER + q_ind.*valDELTIME;
tmp2 = matWVLST.*0;
vcount = zeros(size(matWVLST,1),1);
total_weight = matWVLST(:,1).*0;

for i = 1:3
    dist = sqrt(sum((matWVLST(matWDVE(:,i),:) - matWCENTER).^2,2));
    total_weight(matWDVE(:,i),:) = total_weight(matWDVE(:,i),:) + dist;
end

for i = 1:3
    dist = sqrt(sum((matWVLST(matWDVE(:,i),:) - matWCENTER).^2,2));
    tmp2(matWDVE(:,i),:) = tmp2(matWDVE(:,i),:) + q_ind.*(1 - (dist./total_weight(matWDVE(:,i),:)));
    [gc,gr] = groupcounts(matWDVE(:,i));
    vcount(gr,:) = vcount(gr,:) + gc;
    %     total_weight(matWDVE(:,i),:) = total_weight(matWDVE(:,i),:) + dist;
end
% tmp2 = tmp2./vcount;
% tmp2 = tmp2./(vcount.*total_weight);
% tmp2(isnan(tmp2)) = 0;

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