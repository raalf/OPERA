clear
% clc

%% Preamble

strFILE = 'inputs/2dve.dat';

% [~, strATYPE, vecSYM, ~, ~, ~, valALPHA, valBETA, ~, ~] = fcnOPREAD(strFILE);
strATYPE = 'what'
valALPHA = 0;
valBETA = 0;
vecSYM = [];

matPOINTS(:,:,1) = [0 0 0];
matPOINTS(:,:,3) = [0 1 0];
matPOINTS(:,:,1) = [1 0.5 0];

[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
    matELOC, matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG] = fcnTRIANG(matPOINTS);

vecUINF = fcnUINFWING(valALPHA, 0);

%% Coefficients
% matCOEFF = [ 2 1 0 0 0 0 ];
% matCOEFF = [-45.0959  -41.0083   -5.1624   -0.0404    0.3831    4.1836];
matCOEFF = [-9.1923   45.1363    6.4033    5.5270   -1.4800   -2.9934];

%% Plot

[hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], vecUINF, matROTANG, [3 1 4 4], 'opengl');
% [hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), vecUINF, matROTANG, 'r', 10);
% view([-30 17])

% granularity = 0.1;
% z = -.2:granularity:.2;
% z = -.2:granularity:.2;
% x = -.2:granularity:1.2;

% granularity = 0.01;
% z = -.25:granularity:.25;
% len = length(z);
% y = zeros(len,1) + 0.5;
% x = zeros(len,1) + 0.8;

granularity = 0.3;
z = -2:granularity:2;
len = length(z);
z(end+1) = 0;
y = zeros(len,1) + matCENTER(:,2);% + 0.5;
x = zeros(len,1) + matCENTER(:,1);% + 0.8;
% x = -.2:granularity:1.2;
% y = -.2:granularity:.2;

% granularity = 0.5;
% y = -8:granularity:8;
% z = -1:granularity:1;
% x = -8:granularity:8;

[X,Y,Z] = meshgrid(x,y,z);
fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');

% fpg = [-8 5 1];
% fpg = [0.8 0.5 -0.2];

[s_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);

% q_ind = s_ind + repmat(vecUINF, length(s_ind(:,1)),1);
q_ind = s_ind;

figure(1);
hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3))
hold off







