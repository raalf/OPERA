function [matWELST, matWVLST, matWDVE, valWNELE, matWEIDX, matWPLEX, matWDVECT, matWCENTER, matWROTANG, matWVGRID] = ...
    fcnRELAX5(valDELTIME, valNELE, matCOEFF, matPLEX, valWNELE, matWCOEFF, matWDVE, matWVLST, matWPLEX, valWSIZE, ...
    matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST, matWEIDX, vecDVESYM, vecWDVESYM, vecWSYM, matWVGRID, vecWDVEFLIP, valPRESTEPS)

% Not moving some vertices
% oldest_edge = matWEIDX((1:valWSIZE) + valWSIZE,3); % Trailing edge of oldest wake row
% dont_move = unique([matWELST(vecWLE,:); matWELST(vecWTE,:); matWELST(oldest_edge,:)]);
% dont_move = unique([matWELST(vecWLE,:); matWELST(oldest_edge,:)]);
dont_move = unique([matWELST(vecWLE,:)]);
move = true(size(matWVLST,1),1);
move(dont_move) = false;

presteps = matWVGRID((end-valPRESTEPS+1):end,:);
move(presteps) = false;

%% Getting velocities at wake vertices
tmp2 = zeros(size(matWVLST));
tmp2(move,:) = fcnSDVEVEL(matWVLST(move,:), valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM, [], 1e-3) + fcnSDVEVEL(matWVLST(move,:), valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, vecWDVESYM, [], 1e-3);

tmp2(presteps,:) = tmp2(repmat(matWVGRID(end-valPRESTEPS,:), size(presteps,1), 1),:);

f = matWVGRID(1:end - 2,:);
a = matWVGRID(3:end,:);

ff = [matWVGRID(1,:); f(1:end-1,:)];
aa = [a(2:end,:); matWVGRID(end,:);];

% if size(matWVGRID,1) > 2
%     tmp2(matWVGRID(2:end-1,:),:) = (tmp2(ff,:) + 2.*tmp2(f,:) + 2.*tmp2(matWVGRID(2:end-1,:),:) + 2.*tmp2(a,:) + tmp2(aa,:))./8;
% end

hold on
quiver3(matWVLST(:,1), matWVLST(:,2), matWVLST(:,3), tmp2(:,1), tmp2(:,2), tmp2(:,3));
hold off

%% Moving
% Moving vertices
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


