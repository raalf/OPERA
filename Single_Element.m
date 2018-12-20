clear
% clc

%% Preamble

strFILE = 'inputs/2dve.dat';

% [~, strATYPE, vecSYM, ~, ~, ~, valALPHA, valBETA, ~, ~] = fcnOPREAD(strFILE);
strATYPE = 'what';
valALPHA = 0;
valBETA = 0;
vecSYM = [];

% matPOINTS(:,:,1) = [0 -0.5 0];
% matPOINTS(:,:,2) = [0  0.5 0];
% matPOINTS(:,:,3) = [1  0   0];

matPOINTS(:,:,1) = [0 0 0];
matPOINTS(:,:,2) = [0  0.5 0];
matPOINTS(:,:,3) = [0.5  0.5 0];


[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
    matELOC, matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG, matCONTROL] = fcnTRIANG(matPOINTS);

vecUINF = fcnUINFWING(valALPHA, 0);

%% Coefficients
matCOEFF = [0 0 0 1 0];

%% Plot

% vecUINF = [cosd(10) 0 sind(10)]

[hFig1] = fcnPLOTBODY(1, matDVE, valNELE, matVLST, matELST, matDVECT, matCONTROL, matPLEX, [], vecUINF, matROTANG, [3 1 4 4], 'opengl');
% [hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), vecUINF, matROTANG, 'r', 10);
% view([-30 17])

granularity = 0.05;
y = [-0.2:granularity:0.8];
x = [0.05, matCENTER(:,1), 0.4];
z = [-0.2:granularity:0.2];

% granularity = 0.025;
% x = [-0.2:granularity:0.8];
% y = [0.2, matCENTER(:,2), 0.45];
% z = [-0.2:granularity:0.2];



% granularity = 0.025;
% x = [-0.1:granularity:0.6];
% y = [-0.2:granularity:0.6];
% z = x.*0;

% granularity = 0.01;
% z = -.25:granularity:.25;
% len = length(z);
% y = zeros(len,1) + 0.2;
% x = zeros(len,1) + 0.8;

% granularity = 0.5;
% y = -.5:granularity:1.5;
% z = -1:granularity:1;
% x = -1.5:granularity:1.5;


% granularity = 0.06125;
% y = -.1:granularity:0.6;
% z = -0.5:granularity:0.5;
% x = -0.1:granularity:0.5;

% z(z==0) = [];
[X,Y,Z] = meshgrid(x,y,z);
fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');

% fpg = [matCENTER; matCENTER + [0.1 0 0]];
% fpg = [matCENTER; matCENTER + [0 0 0.1]];
% fpg = matCENTER + [0.1 0.1 0.4];

fpg = [0.2 0 0.02]

% fpg = matCENTER

[q_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCONTROL);

figure(1);
hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3),'b')
hold off







