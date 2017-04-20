clc
clear

fpg = [0 0 1; 0.5 0 0.5];
% fpg = [0 0 1; 0 0 1];
len = length(fpg(:,1));

% matPOINTS(:,:,1) = [0 2.5 0; 0 -2.5, 0];
% matPOINTS(:,:,2) = [0.5 2.5 0; 0 -1.25 0];
% matPOINTS(:,:,3) = [0 1.25 0; 0.5 -2.5 0];

matPOINTS(:,:,1) = [0.5 2.5 0; 0 -2.5, 0];
matPOINTS(:,:,2) = [0 1.25 0; 0 -1.25 0];
matPOINTS(:,:,3) = [0 2.5 0; 0.5 -2.5 0];

% matPOINTS(:,:,1) = [0 1.25 0; 0 -2.5, 0];
% matPOINTS(:,:,2) = [0 2.5 0; 0 -1.25 0];
% matPOINTS(:,:,3) = [0.5 2.5 0; 0.5 -2.5 0];

matCOEFF = [0 1 0 1 0; 0 1 0 1 0];
vecUINF = [0.707 0 0.707];

[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matELOC, matPLEX, matDVECT, matALIGN, matVATT, matVNORM, matCENTER, matROTANG] = fcnTRIANG([], matPOINTS);

[hFig1] = fcnPLOTBODY(1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, vecUINF, matROTANG);

% [q_ind] = fcnINDVEL([1 2], fpg, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, ones(len,1), matROTANG)

[q_ind1] = fcnINDVEL(ones(len,1), fpg, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, ones(len,1), matROTANG);
[q_ind2] = fcnINDVEL(ones(len,1)+1, fpg, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, ones(len,1), matROTANG);



hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind1(:,1), q_ind1(:,2), q_ind1(:,3), 'r')
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind2(:,1), q_ind2(:,2), q_ind2(:,3), 'b')
hold off