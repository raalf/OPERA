function [vecVMU, vecEMU] = fcnVEMU(matVLST, matVATT, matCENTER, matROTANG, matCOEFF, matELST, matEATT, vecTE)
% Gets circulation/Gamma/mu/doublet strength at vertex and edge midpoints
% Without GAPS. Vorticity may not be maintained

% Getting averaged mu values at vertices
vecVMU = nan(size(matVLST,1),1);
for i = 1:size(matVLST,1)
    dves = matVATT(i, ~isnan(matVATT(i,:)))';
    len = length(dves);
    vloc = fcnGLOBSTAR(repmat(matVLST(i,:), len, 1) - matCENTER(dves,:), matROTANG(dves,:));
    circ = sum([0.5.*vloc(:,2).^2 vloc(:,2) 0.5.*vloc(:,1).^2 vloc(:,1) vloc(:,1).*vloc(:,2) ones(size(vloc(:,1)))].*matCOEFF(dves,:),2);
    vecVMU(i) = mean(circ);
end

% Getting averaged mu  values at edge midpoints
vecEMU = nan(size(matELST,1),1);
for i = 1:size(matELST,1)
    dves = matEATT(i, matEATT(i,:) > 0)';
    len = length(dves);
    pt = (matVLST(matELST(i,1),:) + matVLST(matELST(i,2),:))./2;
    vloc = fcnGLOBSTAR(repmat(pt,len,1) - matCENTER(dves,:), matROTANG(dves,:));
    circ = sum([0.5.*vloc(:,2).^2 vloc(:,2) 0.5.*vloc(:,1).^2 vloc(:,1) vloc(:,1).*vloc(:,2) ones(size(vloc(:,1)))].*matCOEFF(dves,:),2);
    vecEMU(i) = mean(circ);
end

% Set to zero leading edge and wing tips
idx = ~all(matEATT,2);
idx(vecTE) = 0;
vecEMU(idx) = 0;
tmp = reshape(matELST(idx,:),[],1,1);
vecVMU(tmp) = 0;

end