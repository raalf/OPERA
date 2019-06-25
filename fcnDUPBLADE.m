function [valNELE, matVLST, matELST, matDVE, matCENTER, matPLEX, vecDVEAREA, ...
    matEATT, matEIDX, vecDVESURFACE, vecDVEFLIP, vecDVESYM, matDVECT, ...
    matROTANG, matVATT, matTEPOINTS, matLEPOINTS] = fcnDUPBLADE(valBLADES, valNELE, matVLST, matELST, matDVE, ...
    matCENTER, matPLEX, vecDVEAREA, matEATT, matEIDX, vecDVESURFACE, vecDVEFLIP, vecDVESYM, matDVECT, matROTANG, matVATT, matTEPOINTS, matLEPOINTS)

theta = (2*pi)/valBLADES;
valNELE = valNELE*valBLADES;

verts = size(matVLST,1);
edges = size(matELST,1);
dves = size(matDVE,1);
tes = size(matTEPOINTS,1);
les = size(matLEPOINTS,1);

for j = 2:valBLADES
    dcmROTORSTEP = angle2dcm(theta*(j-1),0,0,'ZXY');
    
    matVLST = [matVLST; matVLST(1:verts,:)*dcmROTORSTEP];
    matCENTER = [matCENTER; matCENTER(1:dves,:)*dcmROTORSTEP];
    vecDVEAREA = [vecDVEAREA; vecDVEAREA(1:dves)];
    matDVE = [matDVE; matDVE(1:dves,:) + (verts*(j-1))];
    matELST = [matELST; matELST(1:edges,:) + (verts*(j-1))];
    matVATT = [matVATT; matVATT(1:verts,:) + (dves*(j-1)).*(matVATT(1:verts,:) > 0)];
    matEATT = [matEATT; matEATT(1:edges,:) + (dves*(j-1)).*(matEATT(1:edges,:) > 0)];
    matEIDX = [matEIDX; matEIDX(1:dves,:) + (edges*(j-1))];
    vecDVESURFACE = [vecDVESURFACE; vecDVESURFACE(1:dves) + (j-1)];
    vecDVEFLIP = [vecDVEFLIP; vecDVEFLIP(1:dves)];
    
    if ~isempty(vecDVESYM)
        vecDVESYM = [vecDVESYM; vecDVESYM(1:dves,:)];
    end
        
    matTEPOINTS = cat(1, matTEPOINTS, nan(tes,3,2));
    matLEPOINTS = cat(1, matLEPOINTS, nan(les,3,2));
    for i = 1:2
        matTEPOINTS(:,:,i) = [matTEPOINTS(1:(tes*(j-1)),:,i); matTEPOINTS(1:tes,:,i)*dcmROTORSTEP];
        matLEPOINTS(:,:,i) = [matLEPOINTS(1:(les*(j-1)),:,i); matLEPOINTS(1:les,:,i)*dcmROTORSTEP];
    end
    
end

P = permute(reshape(matVLST(matDVE(:,:)',:)', 3, 3, []), [2 1 3]);
DNORM = cross(matVLST(matDVE(:,2),:) - matVLST(matDVE(:,3),:), matVLST(matDVE(:,1),:) - matVLST(matDVE(:,3),:), 2);
DNORM = DNORM./sqrt(sum(DNORM.^2,2));
DNORM(vecDVEFLIP,:) = DNORM(vecDVEFLIP,:).*-1;
[matPLEX, matDVECT, matROTANG] = fcnTRITOLEX(P, DNORM, matCENTER);

end