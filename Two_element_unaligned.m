clear
% clc

%% Preamble

strFILE = 'inputs/2dve.dat';

% [~, strATYPE, vecSYM, ~, ~, ~, valALPHA, valBETA, ~, ~] = fcnOPREAD(strFILE);
strATYPE = 'what';
valALPHA = 10;
valBETA = 0;
vecSYM = [];

matPOINTS(:,:,1) = [0 -0.5 0; 1  -0.5 0];
matPOINTS(:,:,2) = [0  0.5 0; 0  0.5 0];
matPOINTS(:,:,3) = [1 -0.5 0; 1 0.5 0];

[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
    matELOC, matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG, matCONTROL] = fcnTRIANG(matPOINTS);

matUINF = repmat(fcnUINFWING(valALPHA, 0), valNELE, 1);

%% Boundary Conditions
vecTE = 5;
vecLE = 1;
vecTEDVE = 2;
vecSPANDIR = repmat([0 1 0],2,1);
vecSPANDIR = vecSPANDIR - (dot(vecSPANDIR, matDVECT(:,:,3),2)).*matDVECT(:,:,3);
matD = fcnDWING9('2D', matEATT, matPLEX, valNELE, matELOC, matELST, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecSYM, matROTANG, vecSPANDIR, vecTEDVE, matCONTROL);

vecR = 4*pi.*dot(matUINF, matDVECT(:,:,3), 2);
vecR = [zeros(size(matD,1) - 2, 1); vecR];

[matCOEFF] = fcnSOLVED(matD, vecR, valNELE);

% matCOEFF(2,2) = -matCOEFF(2,2)

% matCOEFF(1,:) = [0 0 2.5 1 0];
% matCOEFF(2,:) = [1 -1 0 0 0];

% matCOEFF(1,:) = [0 0 0 0 0];
% matCOEFF(2,:) = [0 0 0 0 0];

%% Plot

[hFig1] = fcnPLOTBODY(1, matDVE, valNELE, matVLST, matELST, matDVECT, matCONTROL, matPLEX, [], matUINF, matROTANG, [3 1 4 4], 'opengl');
[hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), matUINF, matROTANG, 'r', 10);
% view([-30 17])

% granularity = 0.1;
% y = [-3:granularity:.44, .56:granularity:4];
% x = y.*0 + 0.5;
% z = x.*0;

% granularity = 0.1;
% x = [-3:granularity:-0.2, 1.2:granularity:4];
% y = x.*0 + 0.5;
% z = x.*0;

% granularity = 0.1;
% z = -.25:granularity:.25;
% len = length(z);
% y = zeros(len,1) + 0.2;
% x = zeros(len,1) + 0.8;

% granularity = 0.1;
% z = -2:granularity:2;
% len = length(z);
% z(end+1) = 0;
% y = zeros(len,1) + matCENTER(:,2);% + 0.5;
% x = zeros(len,1) + matCENTER(:,1);% + 0.8;
% x = -.2:granularity:1.2;
% y = -.2:granularity:.2;

% granularity = 0.25;
% y = -.5:granularity:3.5;
% z = -1:granularity:1.5;
% x = -2.5:granularity:1.5;

% granularity = 0.25;
% y = -.5:granularity:1.5;
% z = -1:granularity:1;
% x = -2.5:granularity:1.5;
% % z(z==0) = [];

granularity = 0.25;
y = -.625:granularity:.625;
z = -.5:granularity:.5;
x = -.5:granularity:1.5;
z(z==0) = [];

% granularity = 0.25;
% y = -.75:granularity:0.75;
% z = -.5:granularity:.5;
% x = -2.5:granularity:-0.25;
% z(z==0) = [];

% granularity = 1.5;
% y = -4.5:granularity:5.5;
% z = -4.5:granularity:5.5;
% x = -5:granularity:5;

% granularity = 0.1;
% y = 0:granularity:1;
% z = -1:granularity:1;
% x = -0:granularity:1;

% granularity = 0.05;
% y = -0.5:granularity:1.5;
% z = 0:granularity:0; 
% x = -1:granularity:2;

[X,Y,Z] = meshgrid(x,y,z);
fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');
fpg = [fpg; matCENTER];
% fpg = matCENTER + [0 0 0]
% fpg = [0.5 0.5 0]
% fpg = [-1.25 2.25 1.5; -1.0 2.25 1.5];
% fpg = [-1.25 2.25 1.5];

[s_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCONTROL);

q_ind = s_ind + repmat(matUINF(1,:), length(s_ind(:,1)),1);
q_ind = s_ind;

figure(1);
hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3))
% scatter3(fpg(:,1), fpg(:,2), fpg(:,3),500,'xr')
hold off







