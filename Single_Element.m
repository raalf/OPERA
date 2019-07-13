clear
% clc

%% Preamble
% set(groot,'defaultFigureCreateFcn','addToolbarExplorationButtons(gcf)')
% set(groot,'defaultAxesCreateFcn','set(get(gca,''Toolbar''),''Visible'',''off'')')

strFILE = 'inputs/2dve.dat';

% [~, strATYPE, vecSYM, ~, ~, ~, valALPHA, valBETA, ~, ~] = fcnOPREAD(strFILE);
strATYPE = 'what';
valALPHA = 0;
valBETA = 0;
vecSYM = [];

% matPOINTS(:,:,1) = [0 0 0];
% matPOINTS(:,:,2) = [0  0.5 0];
% matPOINTS(:,:,3) = [0.3  0.5 0];

matPOINTS(:,:,1) = [1 0 0];
matPOINTS(:,:,2) = [0  0 0];
% matPOINTS(:,:,3) = [0.9  0.5 0];
xp = 0
% xp = 0.99
% xp = 0.5

matPOINTS(:,:,3) = [xp  -0.5 0];

[TR, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
    matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG, matCONTROL] = fcnTRIANG(matPOINTS,true);

matPOINTS(:,:,3) = [xp  0.5 0];

[TR, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
    matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG, matCONTROL] = fcnTRIANG(matPOINTS, false);

vecDVESYM = false(size(matCENTER, 1), 1);

vecUINF = fcnUINFWING(valALPHA, 0);

matCOEFF = -[0 0 -2 0.17 0 0.056];
% matCOEFF = -[10 0 -20 1 0 1];

%% Plot

[hFig1] = fcnPLOTBODY(1, matDVE, valNELE, matVLST, matELST, matDVECT, matCONTROL, matPLEX, [], vecUINF, matROTANG, [3 1 4 4], 'opengl');
% fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, [], matROTANG, 'r', 10);

% granularity = 0.05;
% x = [-0.1:granularity:0.9];
% y = [-0.1:granularity:0.6];
% z = [-0.1:granularity:0.1];

% % % granularity = 0.00125;
% granularity = 0.04
% x = [-0.1:granularity:0.9];
% y = [-0.1:granularity:0.6];
% % y = [0.01:granularity:0.49];
% z = [-0.15 -0.05 0 0.05 0.15];
% % z = [1 -1];
% z = 0;

% granularity = 0.025
% x = [-0.1:granularity:1];
% y = [-0.2:0.05:0.7];
% z = [-0.05 0 0.05];
% % z = [-0.05 0.05];
% % z = 0;

granularity = 0.05
x = [-0.5:granularity:1.5];
y = [-0.8:granularity:1.3];
z = [-1:granularity:1];
% z = [-0.5 0.5];
z = 0;

% granularity = 0.25
% x = [-1:granularity:2];
% y = [-1:granularity:2];
% z = [-2:granularity:2];
% % z = [1e-1 -1e-1];
% % z = 0;
% % z(z == 0) = []

% granularity = 0.5
% x = [-3:granularity:5];
% y = [-3:granularity:5];
% z = [-4:granularity:5];
% z(z == 0) = []

% granularity = 0.00125;
% granularity = 0.1
% x = [-0.5:granularity:1.5];
% y = [-0.5:granularity:1];
% z = [-0.3:granularity:0.3];

% % granularity = 0.00125;
% granularity = 0.05
% x = [0.05:granularity:0.95];
% % x = 0
% % y = [-1:granularity:1];
% y = [0.05:granularity:0.55];
% % z = [-1:granularity:1];
% z = [1e-5 -1e-5];
% z(z == 0) = []

[X,Y,Z] = meshgrid(x,y,z);
fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');

% num = 100;
% fpg = [linspace(0,xp,num)' linspace(0,0.5,num)' zeros(num,1)];

% fpg = [0.5 0 0; 0.6 0 0; 0.5 0.001 0; 0.6 0.001 0]
% fpg = [0.5 0.2 0; 0.5 0.2 0];

% fpg = [0.3 0.5 0.0; 0.3 0.5 0.0]
% fpg = [0.1 0.1 0.0; 0.1 0.1 0.0]
% fpg = [0.7 0.2 0.0; 0.7 0.2 0.0]
% fpg = [0 1.8 0.8; 0 1.8 0.8]
% fpg = [-0.1 1.3 -0.1; -0.1 1.3 -0.1]
% fpg = [-2 -1 -0.1; -2 -1 -0.1]
% fpg = repmat([1 0 0], 2, 1)

% x1 = (0.2 - 0.5)./(((-0.5)/(1 - xp)));
% fpg = [0 0.2 0; x1 0.2 0];
% vecBOUNDIND = true(size(fpg,1),1);
% vecBOUNDIND = false(size(fpg,1),1);
% fpg = [fpg];

% fpg = [-0.5 0.75 0; -0.5 0.625 0];
% fpg = matCENTER
[q_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCONTROL, vecDVESYM, [], 0);

% vecBOUNDIND = true(size(fpg,1),1);
% [q_ind2] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCONTROL, vecDVESYM, vecBOUNDIND);
% q_ind
% 
% q_ind2 = q_ind + q_ind2;
%%
% figure(1);
% clf(1);
hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3), 'b')
% quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind2(:,1), q_ind2(:,2), q_ind2(:,3), 1, 'm')
hold off

% % figure(2);
% % clf(2);
% plot(z, q_ind(:,3),'-ok')
% grid minor
% box on
% axis tight






