function [matCOEFF] = fcnADJCOEFF(vecVMU, vecEMU, matVLST, matCENTER, matROTANG, matDVE, matCOEFF, matELST, matEIDX, valNELE)

for i = 1:valNELE
    vval = vecVMU(matDVE(i,:));
    tmp1 = fcnDCIRC2(matVLST(matDVE(i,:),:), 1, 1, matROTANG(i,:), matCENTER(i,:));
    
    eloc = (matVLST(matELST(matEIDX(i,:),1),:) + matVLST(matELST(matEIDX(i,:),2),:))./2;
    eval = vecEMU(matEIDX(i,:));
    tmp2 = fcnDCIRC2(eloc, 1, 1, matROTANG(i,:), matCENTER(i,:));
    
    tmpD = [tmp1; tmp2];
    tmpR = [vval; eval];
    
    matCOEFF(i,:) = [tmpD\tmpR]';
end

end