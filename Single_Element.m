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
matPOINTS(:,:,3) = [xp  0.5 0];

[TR, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
    matELOC, matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG, matCONTROL] = fcnTRIANG(matPOINTS, 'SURFACE', []);
vecDVESYM = false(size(matCENTER, 1), 1);

vecUINF = fcnUINFWING(valALPHA, 0);

matCOEFF = -[0 0 -2 0.17 0 0.056];
matCOEFF = -[1 1 1 1 1 1];

%% Plot

[hFig1] = fcnPLOTBODY(1, matDVE, valNELE, matVLST, matELST, matDVECT, matCONTROL, matPLEX, [], vecUINF, matROTANG, [3 1 4 4], 'opengl');
% fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, [], matROTANG, 'r', 10);

granularity = 0.000125;
y = [-0.025:granularity:0.025];
x = [0.2];
z = [-0.025:granularity:0.025].*0;

% % granularity = 0.00125;
% granularity = 0.01
% x = [-0.1:granularity:0.9];
% y = [0.1 0.25 0.4];
% z = [-0.05 0 0.05];
% z = [-0.05 0.05];
% z = 0;

% granularity = 0.00125;
% % granularity = 0.01
% y = [-0.1:granularity:0.9];
% x = [0.1 0.25 0.4];
% % z = [-0.05 0 0.05];
% % z = [-0.05 0.05];
% z = 0;

% granularity = 0.0001;
% x = [0.225:granularity:0.275];
% y = [0.25];
% z = [-0.05:granularity:0.05].*0;

% granularity = 0.05;
% x = [-5:granularity:7];
% y = [-5, matCENTER(:,2), 5];
% z = [-5:granularity:5];

% granularity = 0.05;
% x = 0.2;
% y = [-0.1:granularity:0.1];
% z = 0.0002;

[X,Y,Z] = meshgrid(x,y,z);
fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');

% num = 100;
% fpg = [linspace(0,xp,num)' linspace(0,0.5,num)' zeros(num,1)];

% fpg = [0.5 0 0; 0.6 0 0; 0.5 0.001 0; 0.6 0.001 0]


[q_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCONTROL, vecDVESYM);

%%
figure(1);
hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3), 1, 'b')
hold off






