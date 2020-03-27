function [rs_ie, rs_kk] = fcnRSQUARED(matD, vecR, matCOEFF, boolKINCON, matVLST, matVATT, matCENTER, matROTANG, matELST, matEATT)
% rs_kk = R^2 for kinematic conditions
% rs_ie = R^2 for inter-element boundary conditions

%% ie

% Getting averaged mu values at vertices
vecVMU = nan(size(matVLST,1),1);
for i = 1:size(matVLST,1)
    dves = matVATT(i, ~isnan(matVATT(i,:)))';
    len = length(dves);
    vloc = fcnGLOBSTAR(repmat(matVLST(i,:), len, 1) - matCENTER(dves,:), matROTANG(dves,:));
    circ = sum([0.5.*vloc(:,2).^2 vloc(:,2) 0.5.*vloc(:,1).^2 vloc(:,1) vloc(:,1).*vloc(:,2) ones(size(vloc(:,1)))].*matCOEFF(dves,:),2);
    vecVMU(i) = mean(circ);
    vtx_ss_res(i,1) = sum((circ - vecVMU(i)).^2);
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
    edg_ss_res(i,1) = sum((circ - vecEMU(i)).^2);
end

SS_tot = sum(([vecVMU; vecEMU] - mean([vecVMU; vecEMU])).^2);
SS_res = sum([vtx_ss_res; edg_ss_res]);
rs_ie = 1 - (SS_res/SS_tot);

%% Kin Con
SS_tot = mean((vecR(boolKINCON) - mean(vecR(boolKINCON))).^2);
SS_res = mean((matD(boolKINCON,:)*reshape(matCOEFF',1,[])' - vecR(boolKINCON)).^2);
rs_kk = 1 - (SS_res/SS_tot);

end