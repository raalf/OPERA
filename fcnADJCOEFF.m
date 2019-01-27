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


% [hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], [], matROTANG, [3 1 4 4], 'opengl');
% view([33, 28])
% hold on
% scatter3(matVLST(:,1), matVLST(:,2), -vmu, 'ok');
% scatter3(matCENTER(:,1), matCENTER(:,2), -dmu, '^m');
% pt = (matVLST(matELST(:,1),:) + matVLST(matELST(:,2),:))./2;
% scatter3(pt(:,1), pt(:,2), -emu, 'xr')
% hold off


for i = 1:valNELE
   vloc = fcnGLOBSTAR(matVLST(matDVE(i,:),:) - matCENTER(i,:), repmat(matROTANG(i,:), 3, 1));
   vval = vmu(matDVE(i,:));
   
   eloc = (matVLST(matELST(matEIDX(i,:),1),:) + matVLST(matELST(matEIDX(i,:),2),:))./2;
   eloc = fcnGLOBSTAR(eloc - matCENTER(i,:), repmat(matROTANG(i,:), 3, 1));
   eval = emu(matEIDX(i,:));
   
   dloc = [0 0 0];
   dval = dmu(i);
   
   pts = [vloc; eloc; dloc]; 
   XDATA = pts(:,1); YDATA = pts(:,2);
   ZDATA = [vval; eval; dval]; 
   
   [xData, yData, zData] = prepareSurfaceData( XDATA, YDATA, ZDATA );
   ft = fittype( 'poly22' );
   [fitresult, gof] = fit( [xData, yData], zData, ft);
   newcoeffs = coeffvalues(fitresult);
   order = [6 3 4 2 5 1];
   matCOEFF(i,:) = newcoeffs(order);
   matCOEFF(i,[1,3]) = matCOEFF(i,[1,3]).*2;
end

end