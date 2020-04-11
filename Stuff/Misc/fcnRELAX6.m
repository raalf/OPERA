function [matWELST, matWVLST, matWDVE, valWNELE, matWPLEX, matWDVECT, matWCENTER, matWROTANG] = ...
    fcnRELAX6(flagHVRMOD, valDELTIME, valTIMESTEP, valRPM, valNELE, matCOEFF, matPLEX, valWNELE, matWCOEFF, matWDVE, matWVLST, matWPLEX, valWSIZE, ...
    matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST, matWEIDX, vecDVESYM, vecWDVESYM, vecWSYM, matWVGRID, vecWDVEFLIP, valPRESTEPS, matWE2GRID)
 
presteps = matWVGRID((end - valPRESTEPS+1):end,:);
move = matWVGRID(1:(end - valPRESTEPS),:);

%% Finding induced velocities at all wake control points
q_ind = fcnSDVEVEL(matWCENTER, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM, [], 0) + fcnSDVEVEL(matWCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, vecWDVESYM, [], 0);

tmp2 = zeros(size(matWVLST));
amt = zeros(size(matWVLST,1),1);
for i = 1:valWNELE
    vt = matWDVE(i,:)';
    amt(vt) = amt(vt) + 1;
    tmp2(vt,:) = tmp2(vt,:) + q_ind(i,:);
end
tmp2 = tmp2./amt;

tmp2(presteps,:) = tmp2(repmat(move(end,:), size(presteps,1), 1),:);

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