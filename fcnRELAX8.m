function [matWELST, matWVLST, matWPLEX, matWDVECT, matWCENTER, matWROTANG] = ...
    fcnRELAX8(valTIMESTEP, valDELTIME, valNELE, matCOEFF, matPLEX, valWNELE, matWCOEFF, matWDVE, matWVLST, matWPLEX, valWSIZE, ...
    matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, matWELST, matWEIDX, vecDVESYM, vecWDVESYM, matWVGRID, vecWDVEFLIP, valPRESTEPS, matWDVEGRID)

%% Finding induced velocities at all wake vertices
tmp2 = zeros(size(matWVLST));
sz = 10;
if valTIMESTEP > valPRESTEPS
    for i = 1:sz:(valTIMESTEP-valPRESTEPS)
        if i + sz*2 >(valTIMESTEP-valPRESTEPS)
            idx_rows = [i:(valTIMESTEP-valPRESTEPS)]';
        else
            idx_rows = [i:i+sz]';
        end
        
        q_ind = fcnSDVEVEL(matWCENTER(matWDVEGRID(idx_rows,:),:), valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM, [], 0) + fcnSDVEVEL(matWCENTER(matWDVEGRID(idx_rows,:),:), valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, vecWDVESYM, [], 0);
        
        %     qx = scatteredInterpolant(matWCENTER(matWDVEGRID(idx_rows,:),1), matWCENTER(matWDVEGRID(idx_rows,:),2), matWCENTER(matWDVEGRID(idx_rows,:),3), q_ind(:,1), 'natural', 'nearest');
        %     qy = scatteredInterpolant(matWCENTER(matWDVEGRID(idx_rows,:),1), matWCENTER(matWDVEGRID(idx_rows,:),2), matWCENTER(matWDVEGRID(idx_rows,:),3), q_ind(:,2), 'natural', 'nearest');
        %     qz = scatteredInterpolant(matWCENTER(matWDVEGRID(idx_rows,:),1), matWCENTER(matWDVEGRID(idx_rows,:),2), matWCENTER(matWDVEGRID(idx_rows,:),3), q_ind(:,3), 'natural', 'nearest');
        qx = scatteredInterpolant(matWCENTER(matWDVEGRID(idx_rows,:),1), matWCENTER(matWDVEGRID(idx_rows,:),2), matWCENTER(matWDVEGRID(idx_rows,:),3), q_ind(:,1), 'linear');
        qy = scatteredInterpolant(matWCENTER(matWDVEGRID(idx_rows,:),1), matWCENTER(matWDVEGRID(idx_rows,:),2), matWCENTER(matWDVEGRID(idx_rows,:),3), q_ind(:,2), 'linear');
        qz = scatteredInterpolant(matWCENTER(matWDVEGRID(idx_rows,:),1), matWCENTER(matWDVEGRID(idx_rows,:),2), matWCENTER(matWDVEGRID(idx_rows,:),3), q_ind(:,3), 'linear');
        
        tmp2(matWVGRID(idx_rows,:),:) = [qx(matWVLST(matWVGRID(idx_rows,:),1), matWVLST(matWVGRID(idx_rows,:),2), matWVLST(matWVGRID(idx_rows,:),3)) ...
            qy(matWVLST(matWVGRID(idx_rows,:),1), matWVLST(matWVGRID(idx_rows,:),2), matWVLST(matWVGRID(idx_rows,:),3)) ...
            qz(matWVLST(matWVGRID(idx_rows,:),1), matWVLST(matWVGRID(idx_rows,:),2), matWVLST(matWVGRID(idx_rows,:),3))];
    end
    
    tmp2(matWVGRID, 1) = smoothdata(tmp2(matWVGRID, 1), 1, 'movmedian', 3);
    tmp2(matWVGRID, 2) = smoothdata(tmp2(matWVGRID, 2), 1, 'movmedian', 3);
    tmp2(matWVGRID, 3) = smoothdata(tmp2(matWVGRID, 3), 1, 'movmedian', 3);
    
    tmp2(matWELST(vecWLE,:),:) = tmp2(matWELST(vecWLE,:),:).*0; % Trailing edge of wing
    tmp2(matWVGRID(end,:),:) = tmp2(matWVGRID(end-1,:),:); % Oldest row set equal to second oldest row
end

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