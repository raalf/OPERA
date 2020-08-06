function [valNELE, matVLST, matELST, matDVE, matCENTER, matPLEX, vecDVEAREA, ...
    matEATT, matEIDX, vecDVESURFACE, vecDVEFLIP, matDVECT, ...
    matROTANG, matVATT, matTEPOINTS, matLEPOINTS, matSPANDIR, vecDVEROTOR, vecDVEWING, vecDVESDFLIP] ...
    = fcnREADYROTOR(valROTORS, valNELE, matVLST, matELST, matDVE, ...
    vecDVEAREA, matEATT, matEIDX, vecDVESURFACE, vecDVEFLIP, matVATT, ...
    matTEPOINTS, matLEPOINTS, matSPANDIR, vecROTORBLADES, matROTORHUB, matROTORAXIS, vecDVEWING, vecDVEROTOR, vecDVESDFLIP)

% matELST, matEIDX, matVATT <----------- fix

for n = 1:valROTORS
    idxDVEBLADE = find(vecDVEROTOR == n);
    [idxVLSTBLADE,~,~] = unique(matDVE(idxDVEBLADE,:));
    
    idx_le = ismember(matLEPOINTS(:,:,1), matVLST(idxVLSTBLADE,:), 'rows') & ismember(matLEPOINTS(:,:,2), matVLST(idxVLSTBLADE,:), 'rows');
    idx_te = ismember(matTEPOINTS(:,:,1), matVLST(idxVLSTBLADE,:), 'rows') & ismember(matTEPOINTS(:,:,2), matVLST(idxVLSTBLADE,:), 'rows');
    
    for j = 1:vecROTORBLADES(n)-1
        radBLADE = 2*pi/vecROTORBLADES(n)*j; % radian
        dcmBLADE = angle2dcm(radBLADE,0,0, 'ZXY');
        
        % transformation of rotor from hub plane to xy plane
        dcmXY2HUB = quat2dcm(fcnAXANG2QUAT(vrrotvec([0 0 1], matROTORAXIS(n,:))));
        
        % This is transforming the original single-bladed rotor into hub
        % local, where we will duplicate the blade, then rotate back to the
        % proper location
        POINTS(:,:,1) = ((matVLST(matDVE(idxDVEBLADE,1),:) - matROTORHUB(n,:))/dcmXY2HUB)*dcmBLADE;
        POINTS(:,:,2) = ((matVLST(matDVE(idxDVEBLADE,2),:) - matROTORHUB(n,:))/dcmXY2HUB)*dcmBLADE;
        POINTS(:,:,3) = ((matVLST(matDVE(idxDVEBLADE,3),:) - matROTORHUB(n,:))/dcmXY2HUB)*dcmBLADE;
        
        [~, tmpELST, tmpVLST, tmpDVE, ~, tmpEATT, tmpEIDX, ...
            ~, ~, tmpVATT, ~, ~, ~, ~, ~, ~, ~] = fcnTRIANG(POINTS, vecDVEFLIP(idxDVEBLADE), vecDVESDFLIP(idxDVEBLADE));
        
        vcount = size(matVLST,1);
        dvecount = size(matDVE,1);
        edgecount = size(matELST,1);
        
        % Rotate to rotor axis
        tmpVLST = tmpVLST*dcmXY2HUB + matROTORHUB(n,:);
        
        tmpEATT(tmpEATT ~= 0) = tmpEATT(tmpEATT ~=0) + dvecount;
        matEATT = [matEATT; tmpEATT];
        matELST = [matELST; tmpELST + vcount];
        matEIDX = [matEIDX; tmpEIDX + edgecount];
        matVATT = [matVATT; tmpVATT + dvecount];
        matDVE = [matDVE; tmpDVE + vcount];
        matVLST = [matVLST; tmpVLST];
        
        tmpLEPOINTS(:,:,1) = ((matLEPOINTS(idx_le,:,1) - matROTORHUB(n,:))/dcmXY2HUB)*dcmBLADE;
        tmpLEPOINTS(:,:,1) = tmpLEPOINTS(:,:,1)*dcmXY2HUB + matROTORHUB(n,:);
        
        tmpLEPOINTS(:,:,2) = ((matLEPOINTS(idx_le,:,2) - matROTORHUB(n,:))/dcmXY2HUB)*dcmBLADE;
        tmpLEPOINTS(:,:,2) = tmpLEPOINTS(:,:,2)*dcmXY2HUB + matROTORHUB(n,:);
        
        tmpTEPOINTS(:,:,1) = ((matTEPOINTS(idx_te,:,1) - matROTORHUB(n,:))/dcmXY2HUB)*dcmBLADE;
        tmpTEPOINTS(:,:,1) = tmpTEPOINTS(:,:,1)*dcmXY2HUB + matROTORHUB(n,:);
        
        tmpTEPOINTS(:,:,2) = ((matTEPOINTS(idx_te,:,2) - matROTORHUB(n,:))/dcmXY2HUB)*dcmBLADE;
        tmpTEPOINTS(:,:,2) = tmpTEPOINTS(:,:,2)*dcmXY2HUB + matROTORHUB(n,:);
        
        matLEPOINTS = [matLEPOINTS; tmpLEPOINTS];
        matTEPOINTS = [matTEPOINTS; tmpTEPOINTS];
        
        matSPANDIR = [matSPANDIR; matSPANDIR(idxDVEBLADE,:)*dcmBLADE];
        
        valNELE = valNELE + length(idxDVEBLADE);
        
        vecDVEWING = [vecDVEWING; vecDVEWING(idxDVEBLADE)];
        vecDVEROTOR = [vecDVEROTOR; vecDVEROTOR(idxDVEBLADE)];
        vecDVESURFACE = [vecDVESURFACE; vecDVESURFACE(idxDVEBLADE)];
        vecDVEFLIP = [vecDVEFLIP; vecDVEFLIP(idxDVEBLADE)];
        vecDVESDFLIP = [vecDVESDFLIP; vecDVESDFLIP(idxDVEBLADE)];
        vecDVEAREA = [vecDVEAREA; vecDVEAREA(idxDVEBLADE)];
    end
    
%     % Rotate rotor to axis
%     idxDVEROTOR = vecDVEROTOR == n;
%     idxVLSTROTOR = unique(matDVE(idxDVEROTOR,:));
%     matVLST(idxVLSTROTOR,:)   = matVLST(idxVLSTROTOR,:)  - matROTORHUB(n,:);
%     
%     % Transform rotor from xy plane to hub plane
%     dcmROTORHUB = quat2dcm(fcnAXANG2QUAT(vrrotvec([0 0 1], matROTORAXIS(n,:))));
%     matVLST(idxVLSTROTOR,:)   = matVLST(idxVLSTROTOR,:)  * dcmROTORHUB;
%     matVLST(idxVLSTROTOR,:)   = matVLST(idxVLSTROTOR,:) + matROTORHUB(n,:);
end

matCENTER = (matVLST(matDVE(:,1),:) + matVLST(matDVE(:,2),:) + matVLST(matDVE(:,3),:))./3;
P = permute(reshape(matVLST(matDVE(:,:)',:)', 3, 3, []), [2 1 3]);
DNORM = cross(matVLST(matDVE(:,2),:) - matVLST(matDVE(:,3),:), matVLST(matDVE(:,1),:) - matVLST(matDVE(:,3),:), 2);
DNORM = DNORM./sqrt(sum(DNORM.^2,2));
DNORM(vecDVEFLIP,:) = DNORM(vecDVEFLIP,:).*-1;
[matPLEX, matDVECT, matROTANG] = fcnTRITOLEX(P, DNORM, matCENTER);

end