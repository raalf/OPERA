clear
clc

% profile -memory on

%% Geometry
% matPOINTS = fcnSTLREAD('CAD Geom/naca0015_2d_low.stl');
% matPOINTS = fcnSTLREAD('CAD Geom/naca0015_2d.stl');
% matPOINTS = fcnSTLREAD('CAD Geom/naca0015_2d_high.stl');

% matPOINTS = fcnSTLREAD('CAD Geom/circle_2d_even_low.stl');
% matPOINTS = fcnSTLREAD('CAD Geom/circle_2d_even.stl');
% matPOINTS = fcnSTLREAD('CAD Geom/circle_2d_even_high.stl');

% matPOINTS = fcnSTLREAD('CAD Geom/circle_2d_even_half.stl');
matPOINTS = fcnSTLREAD('CAD Geom/circle_2d_even_half_low.stl');
% matPOINTS = fcnSTLREAD('CAD Geom/circle_2d_even_half_veryhigh.stl');


[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matELOC, matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG] = fcnTRIANG(matPOINTS);

[hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], [], matROTANG, [3 1 4 4], 'opengl');

flagRELAX = 0;
vecSYM = [];

valMAXTIME = 5;
valDELTIME = 0.05;
% valALPHA = atand(1/8)
valALPHA = 0;
matUINF = repmat([1 0 0], valNELE, 1);
% matUINF = repmat([0 0 1], valNELE, 1);

%% D-Matrix Creation
vecTEDVE = [];
vecLEDVE = [];
vecSPANDIR = [];
vecTE = [];
vecLE = [];

matD = fcnDWING9('2D', matEATT, matPLEX, valNELE, matELOC, matELST, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecSYM, matROTANG);

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
% [hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), matUINF, matROTANG, 'r', 20);

hFig1 = figure(1);
hold on
% fpg = matCENTER + (matDVECT(:,:,3)./10000).*-1;
fpg = matCENTER;
q_inds = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);
q_ind = q_inds + matUINF;
fcolor = sqrt(sum(q_ind.^2,2));
fcolor = 1 - fcolor.^(2);

% velocities = [matCOEFF(:,4), matCOEFF(:,2), matCOEFF(:,1).*0];
% fcolor = 1 - (sqrt(sum(velocities.^2,2)) + 1).^2;

p = patch('Faces',matDVE(:,:,1),'Vertices',matVLST,'FaceVertexCData',fcolor,'LineWidth',2);
p.FaceColor = 'flat';
quiver3(fpg(:,1),fpg(:,2),fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3), 'g')
quiver3(fpg(:,1),fpg(:,2),fpg(:,3), q_inds(:,1), q_inds(:,2), q_inds(:,3), 'r')
hold off
colorbar;
grid on
box on
axis equal
view([-42 31])


hFig2 = figure(2);
clf(2)
scatter(fpg(:,1), fcolor, 'xk');
xlabel('Chord Position (m)','FontSize',15);
ylabel('C_p','FontSize',15);
set(gca,'Ydir','reverse')
grid minor
box on
axis tight
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




