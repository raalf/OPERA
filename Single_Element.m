clear
% clc

%% Preamble

strFILE = 'inputs/2dve.dat';

[~, strATYPE, vecSYM, ~, ~, ~, valALPHA, valBETA, ~, ~] = fcnOPREAD(strFILE);

matPOINTS(:,:,1) = [0 0 0];
matPOINTS(:,:,2) = [0 1 0];
matPOINTS(:,:,3) = [1 0.5 0];

[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
    matELOC, matPLEX, matDVECT, matALIGN, matVATT, matVNORM, matCENTER, matROTANG] = fcnTRIANG(matPOINTS);

vecUINF = fcnUINFWING(valALPHA, 0);

%% Coefficients
matCOEFF = [ 0 1 0 0 0 0 ];

%% Plot

[hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], vecUINF, matROTANG, [3 1 4 4], 'opengl');
% [hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), vecUINF, matROTANG, 'r', 10);
view([-30 17])

granularity = 0.1;
y = -0.2:granularity:1.2;
z = -.2:granularity:.2;
x = -.2:granularity:1.2;
granularity = 0.5;
y = -8:granularity:8;
z = -1:granularity:1;
x = -8:granularity:8;



[X,Y,Z] = meshgrid(x,y,z);
fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');

% fpg = [0.100000000000000  -0.200000000000000  -0.200000000000000];

[s_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG, matCENTER);

% q_ind = s_ind + repmat(vecUINF, length(s_ind(:,1)),1);
q_ind = s_ind;

figure(1);
hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3))
hold off







