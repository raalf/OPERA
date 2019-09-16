function [matWELST, matWVLST, matWDVE, valWNELE, matWEIDX, matWPLEX, matWDVECT, matWCENTER, matWROTANG, matWVGRID] = ...
    fcnRELAX7(flagHVRMOD, valDELTIME, valTIMESTEP, valRPM, valNELE, matCOEFF, matPLEX, valWNELE, matWCOEFF, matWDVE, matWVLST, matWPLEX, valWSIZE, ...
    matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST, matWEIDX, vecDVESYM, vecWDVESYM, vecWSYM, matWVGRID, vecWDVEFLIP, valPRESTEPS, matWE2GRID)

presteps = matWVGRID((end - valPRESTEPS+1):end,:);
move = matWVGRID(1:(end - valPRESTEPS),:);
fpg = (matWVLST(move(1:end-1,:),:) + matWVLST(move(2:end,:),:))./2;

%% Getting velocities at wake vertices
tmp2 = zeros(size(matWVLST));
tmp1 = zeros(size(matWELST,1),3);
tmp1(matWE2GRID(1:(end - valPRESTEPS),:),:) = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM, [], 1e-3) + fcnSDVEVEL(fpg, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, vecWDVESYM, [], 1e-3);

if size(move,1) > 2
    tmp2(move(2:end-1,:),:) = (tmp1(matWE2GRID(1:end-1-valPRESTEPS,:),:) + tmp1(matWE2GRID(2:end-valPRESTEPS,:),:))./2;
end
tmp2(move(end,:),:) = tmp1(matWE2GRID(end-valPRESTEPS,:),:);

tmp2(presteps,:) = tmp2(repmat(move(end,:), size(presteps,1), 1),:);

f = matWVGRID(1:end - 2,:);
a = matWVGRID(3:end,:);
if size(matWVGRID,1) > 2
    tmp2(matWVGRID(2:end-1,:),:) = (tmp2(f,:) + 2.*tmp2(matWVGRID(2:end-1,:),:) + tmp2(a,:))./4;
end

% oldest_edge = matWEIDX((1:valWSIZE) + valWSIZE,3); % Trailing edge of oldest wake row
% dont_move = unique([matWELST(vecWLE,:); matWELST(oldest_edge,:)]);
% tmp2(dont_move,:) = tmp2(dont_move,:).*0;

maxRot = 2;
startVel = 6;
if flagHVRMOD && ((valTIMESTEP - valPRESTEPS)*valDELTIME <= (maxRot/(valRPM/60)))
    Vel = -startVel/maxRot*(((valTIMESTEP - valPRESTEPS)*valDELTIME).*(valRPM/60))+startVel;
    tmp2(:,3) = tmp2(:,3) - Vel;
end

dont_move = unique(matWELST(vecWLE,:));
tmp2(dont_move,:) = tmp2(dont_move,:).*0;

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


