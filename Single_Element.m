clear
% clc

%% Preamble

strFILE = 'inputs/2dve.dat';

% [~, strATYPE, vecSYM, ~, ~, ~, valALPHA, valBETA, ~, ~] = fcnOPREAD(strFILE);
strATYPE = 'what';
valALPHA = 0;
valBETA = 0;
vecSYM = [];

matPOINTS(:,:,1) = [0 0 0];
matPOINTS(:,:,2) = [0 1 0];
matPOINTS(:,:,3) = [1 0.5 0];

[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
    matELOC, matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG] = fcnTRIANG(matPOINTS);

vecUINF = fcnUINFWING(valALPHA, 0);

%% Coefficients
matCOEFF = [1 0 0 0 0 0 ];
% matCOEFF = [5.3047    0.0219    5.8048   -0.0268   -0.8324   -5.0658];
% matCOEFF = [-4.0788    1.5006   -2.0544    3.6226   -2.2158   -1.0943];

%% Plot

[hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], vecUINF, matROTANG, [3 1 4 4], 'opengl');
% [hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), vecUINF, matROTANG, 'r', 10);
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
% y = -.5:granularity:1.5;
% z = -1:granularity:1.5;
% x = -2.5:granularity:1.5;

granularity = 0.25;
y = -.5:granularity:1.5;
z = 0.25:granularity:1.5;
x = -2.5:granularity:1.5;

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

% fpg = matCENTER + [0 0 0]

% fpg = [-1 1 0.9];
% fpg = [0.8 0.2 0];
% fpg = [-0.5 -0.25 0]
% fpg = [-2 -1 0]
% fpg = [-0.8 -0.4 0]
% fpg = [0.5 0.95 0]

% fpg = [-2.25 0.5 0.25];
% fpg = [0.5 0.8 0.5];

[s_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);

% q_ind = s_ind + repmat(vecUINF, length(s_ind(:,1)),1);
q_ind = s_ind;

figure(1);
hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3))
hold off







