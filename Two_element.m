clear
% clc

%% Preamble

strFILE = 'inputs/2dve.dat';

% [~, strATYPE, vecSYM, ~, ~, ~, valALPHA, valBETA, ~, ~] = fcnOPREAD(strFILE);
strATYPE = 'what';
valALPHA = 10;
valBETA = 0;
vecSYM = [];

matPOINTS(:,:,1) = [0 -0.5 0; 0  0.5 0];
matPOINTS(:,:,2) = [0  0.5 0; 1  0.5 0];
matPOINTS(:,:,3) = [1 -0.5 0; 1 -0.5 0];

[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
    matELOC, matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG, matCONTROL] = fcnTRIANG(matPOINTS);

matUINF = repmat(fcnUINFWING(valALPHA, 0), valNELE, 1);

%% Coefficients
matCOEFF = zeros(2,6);

%% Kinematic conditions at vertices
% Flow tangency is to be enforced at all control points on the surface HDVEs
% In the D-Matrix, dot (a1,a2,b1,b2,c3) of our influencing HDVE with the normal of the point we are influencing on

% Points we are influencing
fpg = matCENTER;
normals = matDVECT(:,:,3);

% fpg = matCENTER + (matDVECT(:,:,3)./10000).*-1;

% List of DVEs we are influencing from (one for each of the above fieldpoints)
len = length(fpg(:,1));
dvenum = reshape(repmat(1:valNELE,len,1),[],1);
dvetype = ones(size(dvenum));

fpg = repmat(fpg,valNELE,1);

[infl_glob] = fcnHDVEIND(dvenum, dvetype, fpg, matPLEX, matROTANG, matCONTROL);

normals = repmat(normals,valNELE,1); % Repeated so we can dot all at once

% Dotting a1, a2, b1, b2, c3 with the normals of the field points
temp60 = [dot(permute(infl_glob(:,1,:),[3 1 2]),normals,2) dot(permute(infl_glob(:,2,:),[3 1 2]),normals,2) dot(permute(infl_glob(:,3,:),[3 1 2]),normals,2) dot(permute(infl_glob(:,4,:),[3 1 2]),normals,2) dot(permute(infl_glob(:,5,:),[3 1 2]),normals,2) dot(permute(infl_glob(:,6,:),[3 1 2]),normals,2)];

% Reshaping and inserting into the bottom of the D-Matrix
rows = [1:len]';

king_kong = zeros(len, valNELE*6);
king_kong(rows,:) = reshape(permute(reshape(temp60',6,[],valNELE),[2 1 3]),[],6*valNELE,1);
matD = [king_kong(:,4) king_kong(:,8)]

vecR = 4*pi.*dot(matUINF, matDVECT(:,:,3), 2);
coeff = matD\vecR

matCOEFF(1,4) = coeff(1)
matCOEFF(2,2) = coeff(2)


%% Plot

[hFig1] = fcnPLOTBODY(1, matDVE, valNELE, matVLST, matELST, matDVECT, matCONTROL, matPLEX, [], matUINF, matROTANG, [3 1 4 4], 'opengl');
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
% z(z==0) = [];

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
% q_ind = s_ind;

figure(1);
hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3))
% scatter3(fpg(:,1), fpg(:,2), fpg(:,3),500,'xr')
hold off







