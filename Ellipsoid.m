clear
clc

% profile -memory on

%% Geometry

a = 1;
b = 1;
c = 1;

[x,y,z] = ellipsoid(0,0,0,b,a,c,15);
[V,S] = alphavol([x(:), y(:), z(:)]);

TR = triangulation(S.bnd, [x(:), y(:), z(:)]);
matPOINTS = permute(reshape(TR.Points([TR.ConnectivityList(:,1) TR.ConnectivityList(:,3) TR.ConnectivityList(:,2)],:)',3,[],3),[2 1 3]);

[~, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
    matELOC, matPLEX, matDVECT, matALIGN, matVATT, matVNORM, matCENTER, matROTANG] = fcnTRIANG(matPOINTS);

[hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], [], matROTANG, [3 1 4 4], 'opengl');

flagRELAX = 0;
vecSYM = [];

valMAXTIME = 5;
valDELTIME = 0.05;
% valALPHA = atand(1/8)
valALPHA = 0;

%% D-Matrix Creation
vecTEDVE = [];
vecLEDVE = [];
vecSPANDIR = [];
vecTE = [];
vecLE = [];

matD = fcnDWING9([], matEATT, matPLEX, valNELE, matELOC, matELST, matALIGN, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecLEDVE, vecSYM, matVATT, vecTEDVE, vecSPANDIR, matROTANG, matVNORM);

valDLEN = length(matD);

%% Preparing to timestep
valTIMESTEP = 0;
matWADJE = [];
matWELST = [];
matWDVE = [];
valWNELE = [];
matWEATT = [];
matWEIDX = [];
matWELOC = [];
matWPLEX = [];
matWDVECT = [];
matWALIGN = [];
matWVATT = [];
matWVNORM = [];
matWCENTER = [];
matWAKEGEOM = [];
valWSIZE = [];
matWCOEFF = [];
matWVLST = [];

vecUINF = fcnUINFWING(valALPHA, 0);
matUINF = repmat(vecUINF,valNELE,1);

% Building wing resultant
vecR = fcnRWING([], valDLEN, 0, matELST, matCENTER, matDVECT, matUINF(1,:), vecLE, vecLEDVE, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, [], matVNORM, matVLST);

% Solving for wing coefficients
[matCOEFF] = fcnSOLVED(matD, vecR, valNELE);

%% Plot
% [hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], matUINF, matROTANG, [3 1 4 4], 'opengl');
% [hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), matUINF, matROTANG, 'r', 5);

% q_inds = fcnSDVEVEL(matVLST, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG, matCENTER);
% q_ind = q_inds + matUINF(1,:);
% matVLSTCP = sqrt(sum(q_ind.^2,2));
% matVLSTCP = 1 - sqrt(sum(q_ind.^2,2)).^2;

hFig1 = figure(1);
hold on

q_inds = fcnSDVEVEL(matCENTER, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG, matCENTER);
q_ind = q_inds + matUINF(1,:);
fcolor = sqrt(sum(q_ind.^2,2));
p = patch('Faces',matDVE(:,:,1),'Vertices',matVLST,'FaceVertexCData',fcolor,'LineWidth',2);
p.FaceColor = 'flat';

quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3), 'g')
hold off

colorbar;
grid on
box on
axis equal



granularity = .25;
x = -3:granularity:3;
% y = -1.2:granularity:1.2;
y = ones(size(x)) - 1;
z = -3:granularity:3;

% granularity = 1;
% x = -30:granularity:30;
% % y = -1.2:granularity:1.2;
% y = ones(size(x)) - 1;
% z = -30:granularity:30;

[X,Y,Z] = meshgrid(x,y,z);
fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');

[s_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG, matCENTER);

q_ind = s_ind + repmat(vecUINF, length(s_ind(:,1)),1);
% q_ind = s_ind;
hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3))
hold off
axis tight




