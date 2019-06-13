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
xp = 0.99
% xp = 0
% xp = 0.5
matPOINTS(:,:,3) = [xp  0.5 0];

[TR, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
    matELOC, matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG, matCONTROL] = fcnTRIANG(matPOINTS, 'SURFACE', [], false);
vecDVESYM = false(size(matCENTER, 1), 1);

vecUINF = fcnUINFWING(valALPHA, 0);

matCOEFF = -[0 0 -2 0.17 0 0.056];
matCOEFF = -[1 1 1 1 1 1];

%% Plot

[hFig1] = fcnPLOTBODY(1, matDVE, valNELE, matVLST, matELST, matDVECT, matCONTROL, matPLEX, [], vecUINF, matROTANG, [3 1 4 4], 'opengl');
% fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, [], matROTANG, 'r', 10);

granularity = 0.00125;
y = [0.5];
x = [0.99];
z = [-0.25:granularity:0.25];

% % granularity = 0.00125;
% granularity = 0.1
% x = [-0.1:granularity:0.9];
% y = [-0.1:granularity:0.6];
% % y = [0.01:granularity:0.49];
% z = [-0.05 0 0.05];
% z = [-0.05 0.05];
% % z = 0;

% % granularity = 0.00125;
% granularity = 0.1
% % x = [-0.1:granularity:0.9];
% x = [0.6:granularity:0.9];
% 
% granularity = 0.1
% x = [0.6:granularity:0.9];
% y = [0.2];
% % y = [0.01:granularity:0.49];
% z = [-0.05 0 0.05];
% % z = [-0.05 0.05];
% % z = 0;

% granularity = 0.00125;
granularity = 0.0625
x = [-1:granularity:2];
% x = 0
% y = [-1:granularity:1];
y = [-1:granularity:2];
z = [-0.3:granularity:0.3];
% z = [1e-5 -1e-5];
% z = 0;
% z(z == 0) = []

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


[q_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCONTROL, vecDVESYM);



%%
% figure(1);
% clf(1);
hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3), 1, 'b')
hold off

% % figure(2);
% % clf(2);
% plot(z, q_ind(:,3),'-ok')
% grid minor
% box on
% axis tight






