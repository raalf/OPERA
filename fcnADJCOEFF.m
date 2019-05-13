function matCOEFF = fcnADJCOEFF(matVLST, matVATT, matCENTER, matROTANG, matDVE, matCOEFF, matELST, matEATT, matEIDX, vecTE, valNELE, matERROR)

% Getting averaged mu values at vertices
vmu = nan(size(matVLST,1),1);
for i = 1:size(matVLST,1)
    dves = matVATT(i, ~isnan(matVATT(i,:)))';
    len = length(dves);
    vloc = fcnGLOBSTAR(repmat(matVLST(i,:), len, 1) - matCENTER(dves,:), matROTANG(dves,:));
    circ = sum([0.5.*vloc(:,2).^2 vloc(:,2) 0.5.*vloc(:,1).^2 vloc(:,1) vloc(:,1).*vloc(:,2) ones(size(vloc(:,1)))].*matCOEFF(dves,:),2);
    vmu(i) = mean(circ);
end

% Getting averaged mu  values at edge midpoints
emu = nan(size(matELST,1),1);
for i = 1:size(matELST,1)
    dves = matEATT(i, matEATT(i,:) > 0)';
    len = length(dves);
    pt = (matVLST(matELST(i,1),:) + matVLST(matELST(i,2),:))./2;
    vloc = fcnGLOBSTAR(repmat(pt,len,1) - matCENTER(dves,:), matROTANG(dves,:));
    circ = sum([0.5.*vloc(:,2).^2 vloc(:,2) 0.5.*vloc(:,1).^2 vloc(:,1) vloc(:,1).*vloc(:,2) ones(size(vloc(:,1)))].*matCOEFF(dves,:),2);
    emu(i) = mean(circ);
end

% Set to zero leading edge and wing tips
idx = ~all(matEATT,2);
idx(vecTE) = 0;
emu(idx) = 0;
tmp = reshape(matELST(idx,:),[],1,1);
vmu(tmp) = 0;

% Getting mu  values at DVE control points
dmu = matCOEFF(:,6);

for i = 1:valNELE
    vval = vmu(matDVE(i,:));
    tmp1 = fcnDCIRC2(matVLST(matDVE(i,:),:), 1, 1, matROTANG(i,:), matCENTER(i,:));
        
    eloc = (matVLST(matELST(matEIDX(i,:),1),:) + matVLST(matELST(matEIDX(i,:),2),:))./2;
    eval = emu(matEIDX(i,:));
    tmp2 = fcnDCIRC2(eloc, 1, 1, matROTANG(i,:), matCENTER(i,:));
    
    dloc = matCENTER(i,:);
    dval = dmu(i);
    tmp3 = fcnDCIRC2(dloc, 1, 1, matROTANG(i,:), matCENTER(i,:));
    
%     tmpD = [tmp1; tmp2; tmp3];
%     tmpR = [vval; eval; dval];
    
    tmpD = [tmp1; tmp2];
    tmpR = [vval; eval];

    matCOEFF(i,:) = [tmpD\tmpR]';
end

end