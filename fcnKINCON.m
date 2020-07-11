function [matKINCON_P, vecKINCON_DVE] = fcnKINCON(valNELE, matDVE, matCENTER, matVLST, vecDVEAREA)

tmp = find(vecDVEAREA > .5e-4);
vecKINCON_DVE = [1:valNELE repmat(tmp',1,3)]';
dist = 0.8;
matKINCON_P = [matCENTER; ...
    dist.*(matVLST(matDVE(tmp,1),:) - matCENTER(tmp,:)) + matCENTER(tmp,:); ...
    dist.*(matVLST(matDVE(tmp,2),:) - matCENTER(tmp,:)) + matCENTER(tmp,:); ...
    dist.*(matVLST(matDVE(tmp,3),:) - matCENTER(tmp,:)) + matCENTER(tmp,:)];

end