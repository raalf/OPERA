clear
clc

% profile -memory on

%% Geometry

matPOINTS = fcnSTLREAD('CAD Geom/naca0015_low.stl');

[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matELOC, matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG] = fcnTRIANG(matPOINTS);

[hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], [], matROTANG, [3 1 4 4], 'opengl');

flagRELAX = 0;
vecSYM = [];

valMAXTIME = 5;
valDELTIME = 0.05;
% valALPHA = atand(1/8)
valALPHA = 0;
matUINF = repmat([1 0 0], valNELE, 1);

%% D-Matrix Creation
vecTEDVE = [];
vecLEDVE = [];
vecSPANDIR = [];
vecTE = [];
vecLE = [];

matD = fcnDWING9([], matEATT, matPLEX, valNELE, matELOC, matELST, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecSYM, matROTANG);

valDLEN = length(matD);

%% Preparing to timestep
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
matWROTANG = [];

% Building wing resultant
vecR = fcnRWING(valDLEN, 0, matCENTER, matDVECT, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER);

% Solving for wing coefficients
[matCOEFF] = fcnSOLVED(matD, vecR, valNELE);


%% Plot
[hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], matUINF, matROTANG, [3 1 4 4], 'opengl');
[hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), matUINF, matROTANG, 'r', 20);

% %% Plot
% 
% hFig1 = figure(1);
% hold on
% q_inds = fcnSDVEVEL(matCENTER, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);
% q_ind = q_inds + matUINF(1,:);
% fcolor = sqrt(sum(q_ind.^2,2));
% p = patch('Faces',matDVE(:,:,1),'Vertices',matVLST,'FaceVertexCData',fcolor,'LineWidth',2);
% p.FaceColor = 'flat';
% % quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3), 'g')
% hold off
% colorbar;
% grid on
% box on
% axis equal
% 
% granularity = .5;
% x = -1.5:granularity:1.5;
% % y = -3:granularity:3;
% y = -1.5:granularity:1.5;
% z = -1.5:granularity:1.5;
% [X,Y,Z] = meshgrid(x,y,z);
% 
% s_ind = fcnSDVEVEL([X(:) Y(:) Z(:)], valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);
% q_ind = s_ind + repmat(matUINF(1,:), length(s_ind(:,1)),1);
% Xq = reshape(q_ind(:,1), size(X));
% Yq = reshape(q_ind(:,2), size(Y));
% Zq = reshape(q_ind(:,3), size(Z));
% 
% [Xs,Ys,Zs] = meshgrid(-1.5,y,z);
% hold on
% streamline(X,Y,Z,Xq,Yq,Zq,Xs,Ys,Zs);
% % quiver3(X(:),Y(:),Z(:),q_ind(:,1),q_ind(:,2),q_ind(:,3));
% hold off




