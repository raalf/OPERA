clc
clear

matPOINTS(:,:,1) = [0 0 0];
matPOINTS(:,:,2) = [0 2 0];
matPOINTS(:,:,3) = [1 0 0];

matCOEFF = [0 0 0 1 0 0];

fpg = [10 20 30];

% [TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
%     matELOC, matPLEX, matDVECT, matALIGN, matVATT, matVNORM, matCENTER, matROTANG] = fcnTRIANG(matPOINTS);

xi_1 = matPOINTS(1,1,1);
xi_2 = matPOINTS(1,1,2);
xi_3 = matPOINTS(1,1,3);
xi_12 = xi_1;

eta_1 = matPOINTS(1,2,1);
eta_2 = matPOINTS(1,2,2);
eta_3 = matPOINTS(1,2,3);

xi_p = fpg(1);
eta_p = fpg(2);
zeta_p = fpg(3);

A1 = matCOEFF(1);
A2 = matCOEFF(2);
B1 = matCOEFF(3);
B2 = matCOEFF(4);
C2 = matCOEFF(5);
C3 = matCOEFF(6);



w_eta = (A1*zeta_p*sign(xi_3 - xi_p)*sign(eta_1 - eta_p - (eta_1*xi_1)/(xi_1 - xi_3) + (eta_1*xi_3)/(xi_1 - xi_3) + (eta_3*xi_1)/(xi_1 - xi_3) - (eta_3*xi_3)/(xi_1 - xi_3))*log((eta_3*xi_1 - eta_3*xi_3 - eta_p*xi_1 + eta_p*xi_3)/(xi_1 - xi_3)))/(eta_1/(xi_1 - xi_3) - eta_3/(xi_1 - xi_3))
