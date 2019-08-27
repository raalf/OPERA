function [matHINGELOC, matHINGEDIR, matVLST, matCENTER, matPLEX, matDVECT, matROTANG] ...
    = fcnSETHANG(matHINGEANG, matHINGELOC, matHINGEDIR, matVLST, matDVE, vecDVESURFACE, valBLADES, vecDVEFLIP)

for j = 1:valBLADES
    
    % FLAP
    % Shift pitch and hunt axis and locations, and all dves on this blade
    tmpR = fcnROTMAX(matHINGEDIR(j,:,1), -matHINGEANG(j,1));
    matHINGEDIR(j,:,2) = matHINGEDIR(j,:,2)*tmpR;
    matHINGEDIR(j,:,3) = matHINGEDIR(j,:,3)*tmpR;
    
    matHINGELOC(j,:,2) = (matHINGELOC(j,:,2) - matHINGELOC(j,:,1))*tmpR + matHINGELOC(j,:,1);
    matHINGELOC(j,:,3) = (matHINGELOC(j,:,3) - matHINGELOC(j,:,1))*tmpR + matHINGELOC(j,:,1);   
    
    matVLST(matDVE(vecDVESURFACE == j,:),:) = (matVLST(matDVE(vecDVESURFACE == j,:),:) - matHINGELOC(j,:,1))*tmpR + matHINGELOC(j,:,1);
    
    % PITCH
    tmpR = fcnROTMAX(matHINGEDIR(j,:,2), -matHINGEANG(j,2));
    matHINGEDIR(j,:,3) = matHINGEDIR(j,:,3)*tmpR;
    matHINGELOC(j,:,3) = (matHINGELOC(j,:,3) - matHINGELOC(j,:,2))*tmpR + matHINGELOC(j,:,2);   
    
    matVLST(matDVE(vecDVESURFACE == j,:),:) = (matVLST(matDVE(vecDVESURFACE == j,:),:) - matHINGELOC(j,:,2))*tmpR + matHINGELOC(j,:,2);    
    
    % HUNT   
    matVLST(matDVE(vecDVESURFACE == j,:),:) = (matVLST(matDVE(vecDVESURFACE == j,:),:) - matHINGELOC(j,:,3))*tmpR + matHINGELOC(j,:,3);        
    
end

matCENTER = (matVLST(matDVE(:,1),:) + matVLST(matDVE(:,2),:) + matVLST(matDVE(:,3),:))./3;
P = permute(reshape(matVLST(matDVE(:,:)',:)', 3, 3, []), [2 1 3]);
DNORM = cross(matVLST(matDVE(:,2),:) - matVLST(matDVE(:,3),:), matVLST(matDVE(:,1),:) - matVLST(matDVE(:,3),:), 2);
DNORM = DNORM./sqrt(sum(DNORM.^2,2));
DNORM(vecDVEFLIP,:) = DNORM(vecDVEFLIP,:).*-1;
[matPLEX, matDVECT, matROTANG] = fcnTRITOLEX(P, DNORM, matCENTER);


% hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], [], matROTANG, [], 'opengl');
% hold on
% quiver3(matHINGELOC(:,1,1), matHINGELOC(:,2,1), matHINGELOC(:,3,1), matHINGEDIR(:,1,1), matHINGEDIR(:,2,1), matHINGEDIR(:,3,1), 'b');
% quiver3(matHINGELOC(:,1,2), matHINGELOC(:,2,2), matHINGELOC(:,3,2), matHINGEDIR(:,1,2), matHINGEDIR(:,2,2), matHINGEDIR(:,3,2), 'r');
% quiver3(matHINGELOC(:,1,3), matHINGELOC(:,2,3), matHINGELOC(:,3,3), matHINGEDIR(:,1,3), matHINGEDIR(:,2,3), matHINGEDIR(:,3,3), 'm');
% hold off

end