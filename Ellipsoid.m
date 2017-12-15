clear
clc

%% Geometry

a = 2;
b = 1.5;
c = 1;

[x,y,z] = ellipsoid(0,0,0,b,a,c,20);
[V,S] = alphavol([x(:), y(:), z(:)]);

TR = triangulation(S.bnd, [x(:), y(:), z(:)]);
matPOINTS = permute(reshape(TR.Points([TR.ConnectivityList(:,1) TR.ConnectivityList(:,3) TR.ConnectivityList(:,2)],:)',3,[],3),[2 1 3]);

[~, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
    matELOC, matPLEX, matDVECT, matALIGN, matVATT, matVNORM, matCENTER, matROTANG, matVSCOMB] = fcnTRIANG(matPOINTS);

flagRELAX = 0;
vecSYM = [];

valMAXTIME = 5;
valDELTIME = 0.05;
valALPHA = 0;

%% D-Matrix Creation
vecTEDVE = [];
vecLEDVE = [];
vecSPANDIR = [];
vecTE = [];
vecLE = [];

matD = fcnDWING9([], matEATT, matPLEX, valNELE, matELOC, matELST, matALIGN, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecLEDVE, vecSYM, matVATT, vecTEDVE, vecSPANDIR, matROTANG, matVNORM, matVSCOMB);

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
vecR = fcnRWING([], valDLEN, 0, matELST, matCENTER, matDVECT, matUINF, vecLE, vecLEDVE, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, [], matVNORM, matVLST);

% Solving for wing coefficients
[matCOEFF] = fcnSOLVED(matD, vecR, valNELE);

%% Plot

[hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], matUINF, matROTANG, [3 1 4 4], 'opengl');
% [hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), matUINF, matROTANG, 'r', 5);

% hold on
% q_inds = fcnSDVEVEL(matCENTER, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG, [], matCENTER);
% q_ind = q_inds + matUINF;
% quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3), 0.25, 'g')
% hold off


granularity = .5;
x = -3:granularity:3;
% y = -1.2:granularity:1.2;
y = ones(size(x)) - 1;
z = -3:granularity:3;
[X,Y,Z] = meshgrid(x,y,z);
fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');

[s_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG, matVSCOMB, matCENTER);

q_ind = s_ind + repmat(vecUINF, length(s_ind(:,1)),1);
% q_ind = s_ind;
hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3))
hold off
axis tight

granularity = 1;
x = -3:granularity:3;
y = -3:granularity:3;
z = -3:granularity:3;
[X,Y,Z] = meshgrid(x,y,z);

s_ind = fcnSDVEVEL([X(:) Y(:) Z(:)],valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG, matVSCOMB, matCENTER);
q_ind = s_ind + repmat(vecUINF, length(s_ind(:,1)),1);
Xq = reshape(q_ind(:,1), size(X));
Yq = reshape(q_ind(:,2), size(Y));
Zq = reshape(q_ind(:,3), size(Z));

[Xs,Ys,Zs] = meshgrid(-3,y,z);
hold on
streamline(X,Y,Z,Xq,Yq,Zq,Xs,Ys,Zs);
% quiver3(X(:),Y(:),Z(:),q_ind(:,1),q_ind(:,2),q_ind(:,3));
hold off






