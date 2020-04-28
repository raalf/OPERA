function [matWELST, matWVLST, matWDVE, valWNELE, matWEIDX, matWPLEX, matWDVECT, matWCENTER, matWROTANG] = ...
    fcnRELAX6(valTIMESTEP, valDELTIME, valNELE, matCOEFF, matPLEX, valWNELE, matWCOEFF, matWDVE, matWVLST, matWPLEX, valWSIZE, ...
    matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, matWELST, matWEIDX, vecDVESYM, vecWDVESYM, matWVGRID, vecWDVEFLIP, valPRESTEPS)

dvegrid = flipud(repmat([1:valWSIZE, valWSIZE*2], valTIMESTEP, 1) + [0:(valWSIZE*2):(valWSIZE*valTIMESTEP*2 - 1)]');

%% Finding induced velocities at all wake vertices
q_ind = fcnSDVEVEL(matWCENTER, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM, [], 0) + fcnSDVEVEL(matWCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, vecWDVESYM, [], 0);

qx = scatteredInterpolant(matWCENTER(:,1), matWCENTER(:,2), matWCENTER(:,3), q_ind(:,1), 'natural', 'nearest');
qy = scatteredInterpolant(matWCENTER(:,1), matWCENTER(:,2), matWCENTER(:,3), q_ind(:,2), 'natural', 'nearest');
qz = scatteredInterpolant(matWCENTER(:,1), matWCENTER(:,2), matWCENTER(:,3), q_ind(:,3), 'natural', 'nearest');

tmp2 = [qx(matWVLST(:,1), matWVLST(:,2), matWVLST(:,3)) qy(matWVLST(:,1), matWVLST(:,2), matWVLST(:,3)) qz(matWVLST(:,1), matWVLST(:,2), matWVLST(:,3))];

tmp2(matWVGRID, 1) = smoothdata(tmp2(matWVGRID, 1), 1);
tmp2(matWVGRID, 2) = smoothdata(tmp2(matWVGRID, 2), 1);
tmp2(matWVGRID, 3) = smoothdata(tmp2(matWVGRID, 3), 1);

tmp2(matWELST(vecWLE,:),:) = tmp2(matWELST(vecWLE,:),:).*0; % Trailing edge of wing

% % tmp2(matWELST(vecWTE,:),:) = tmp2(matWELST(vecWTE,:),:).*0; % Trailing edge of LE wake element row
% oldest_edge = matWEIDX((1:valWSIZE) + valWSIZE,3); % Trailing edge of oldest wake row
% tmp2(matWELST(oldest_edge,:),:) = tmp2(matWELST(oldest_edge,:),:).*0;
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