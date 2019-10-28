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

matPOINTS(:,:,3) = [xp 0.5 0];

[TR, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
    matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG, matCONTROL] = fcnTRIANG(matPOINTS, false);

% matPOINTS(:,:,3) = [xp  -0.5 0];
% 
% [TR, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
%     matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG, matCONTROL] = fcnTRIANG(matPOINTS,true);

vecDVESYM = false(size(matCENTER, 1), 1);

vecUINF = fcnUINFWING(valALPHA, 0);

% matCOEFF = [0 0 0   -0.3236    0.0000   -0.1021];
% matCOEFF = -[10 0 -20 1 0 1];
matCOEFF = -[1 1 1 1 1 1];

%% Plot

[hFig1] = fcnPLOTBODY(1, matDVE, valNELE, matVLST, matELST, matDVECT, matCONTROL, matPLEX, [], vecUINF, matROTANG, [3 1 4 4], 'opengl');
% fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, [], matROTANG, 'r', 10);

% granularity = 0.05;
% x = [-0.1:granularity:0.9];
% y = [-0.1:granularity:0.6];
% z = [-0.1:granularity:0.1];

granularity = 0.125;
x = [-0.5:granularity:1.5];
y = [-0.8:granularity:1.3];
% z = [-0.5:granularity.*2:0.5];
z = [-0.5 -0.25 0.25 0.5];
% z = 0;

% granularity = 0.0125
% y = [-0.5:granularity:1.5];
% x = matCENTER(1);
% z = matCENTER(3);

% granularity = 0.000125
% x = [-0.1:granularity:0.1];
% y = 0.25;
% z = 0;

% granularity = 0.5
% x = [-3:granularity:5];
% y = [-3:granularity:5];
% z = [-4:granularity:5];
% z(z == 0) = []
% % z = 0

% granularity = 0.01;
% x = matCENTER(1);
% y = matCENTER(2);
% x = 0;
% y = 0.5;
% x = 0.25;
% y = 0.25;
% x = 1;
% y = 2;
% z = [-1:granularity:1];

[X,Y,Z] = meshgrid(x,y,z);
fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');

% fpg = [0 0.25 0]
% fpg = [0 0 0; 1 0 0; 0.5 0.5 0]
% fpg = [x 0 z];
% fpg = [1 4.75 0; 1 4.5 0]
fpg = [-0.125 -0.675 -0.25; -0.25 -0.675 -0.25]
[q_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCONTROL, vecDVESYM, [], 0);
[q_ind2] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCONTROL, vecDVESYM, [], 0.001);

%%
% figure(2);
% clf(2);
% plot(q_ind(:,3), y, '-ok')
% hold on
% plot(q_ind2(:,3), y, '--b^')
% hold off
% grid minor
% box on
% axis tight

%%
% figure(4);
% clf(4);
% subplot(3,1,1)
% plot(q_ind(:,1), z, '-ok')
% hold on
% plot(q_ind2(:,1), z, '--b^')
% hold off
% grid minor
% box on
% axis tight
% 
% subplot(3,1,2)
% plot(q_ind(:,2), z, '-ok')
% hold on
% plot(q_ind2(:,2), z, '--b^')
% hold off
% grid minor
% box on
% axis tight
% 
% subplot(3,1,3)
% plot(q_ind(:,3), z, '-ok')
% hold on
% plot(q_ind2(:,3), z, '--b^')
% hold off
% grid minor
% box on
% axis tight
% 
% points = fcnGLOBSTAR([x y 0] - matCENTER, matROTANG);
% vort = [matCOEFF(1,3).*points(:,1) + matCOEFF(1,4), matCOEFF(1,1).*points(:,2) + matCOEFF(1,2), points(:,2).*0];

% vecBOUNDIND = true(size(fpg,1),1);
% [q_ind2] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCONTROL, vecDVESYM, vecBOUNDIND);

%%
figure(1);
% clf(1);
hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3), 0,'b')
% quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind2(:,1), q_ind2(:,2), q_ind2(:,3), 1, 'm')
hold off







