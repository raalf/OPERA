function [matWELST, matWVLST, matWDVE, valWNELE, matWEIDX, matWPLEX, matWDVECT, matWCENTER, matWROTANG, matWVGRID] = ...
    fcnRELAX7(valDELTIME, valNELE, matCOEFF, matPLEX, valWNELE, matWCOEFF, matWDVE, matWVLST, matWPLEX, valWSIZE, ...
    matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST, matWEIDX, vecDVESYM, vecWDVESYM, vecWSYM, matWVGRID, vecWDVEFLIP, valPRESTEPS, matWE2GRID)

% Not moving some vertices
% oldest_edge = matWEIDX((1:valWSIZE) + valWSIZE,3); % Trailing edge of oldest wake row
% oldest_verts = unique([matWELST(oldest_edge,:)]);

presteps = matWVGRID((end-valPRESTEPS+1):end,:);

move = matWVGRID(1:(end - valPRESTEPS),:);
fpg = (matWVLST(move(1:end-1,:),:) + matWVLST(move(2:end,:),:))./2;

%% Getting velocities at wake vertices
tmp2 = zeros(size(matWVLST));
tmp1 = zeros(size(matWELST,1),3);
tmp1(matWE2GRID(1:(end - valPRESTEPS),:),:) = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM, [], 1e-14) + fcnSDVEVEL(fpg, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, vecWDVESYM, [], 1e-14);

% hold on
% quiver3(fpg(:,1), fpg(:,2), fpg(:,3), tmp1(matWE2GRID(1:(end - valPRESTEPS),:),1), tmp1(matWE2GRID(1:(end - valPRESTEPS),:),2), tmp1(matWE2GRID(1:(end - valPRESTEPS),:),3),'m');
% hold off

if size(move,1) > 2
    tmp2(move(2:end-1,:),:) = (tmp1(matWE2GRID(1:end-1-valPRESTEPS,:),:) + tmp1(matWE2GRID(2:end-valPRESTEPS,:),:))./2; 
end
tmp2(move(end,:),:) = tmp1(matWE2GRID(end-valPRESTEPS,:),:);

% hold on
% quiver3(matWVLST(:,1), matWVLST(:,2), matWVLST(:,3), tmp2(:,1), tmp2(:,2), tmp2(:,3),'r');
% hold off

% hold on
% quiver3([fpg(:,1); matWVLST(move,1)], [fpg(:,2); matWVLST(move,2)], [fpg(:,3); matWVLST(move,3)], [tmp1(matWE2GRID(1:(end - valPRESTEPS),:),1); tmp2(move,1)], [tmp1(matWE2GRID(1:(end - valPRESTEPS),:),2); tmp2(move,2)], [tmp1(matWE2GRID(1:(end - valPRESTEPS),:),3); tmp2(move,3)]);
% hold off

tmp2(presteps,:) = tmp2(repmat(move(end,:), size(presteps,1), 1),:);

f = matWVGRID(1:end - 2,:);
a = matWVGRID(3:end,:);

ff = [matWVGRID(1,:); f(1:end-1,:)];
aa = [a(2:end,:); matWVGRID(end,:);];

if size(matWVGRID,1) > 2
    tmp2(matWVGRID(2:end-1,:),:) = (tmp2(f,:) + 2.*tmp2(matWVGRID(2:end-1,:),:) + tmp2(a,:))./4;
end

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


