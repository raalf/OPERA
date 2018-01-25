clear
clc

% profile -memory on

%% Geometry

a = 1;
b = 1;
c = 1;

[x,y,z] = ellipsoid(0,0,0,b,a,c,14);
[V,S] = alphavol([x(:), y(:), z(:)]);

TR = triangulation(S.bnd, [x(:), y(:), z(:)]);
matPOINTS = permute(reshape(TR.Points([TR.ConnectivityList(:,1) TR.ConnectivityList(:,3) TR.ConnectivityList(:,2)],:)',3,[],3),[2 1 3]);

[~, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
    matELOC, matPLEX, matDVECT, matALIGN, matVATT, matVNORM, matCENTER, matROTANG] = fcnTRIANG(matPOINTS);

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
vecR = fcnRWING([], valDLEN, 0, matELST, matCENTER, matDVECT, matUINF, vecLE, vecLEDVE, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, [], matVNORM, matVLST);

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
clf(1);
q_inds = fcnSDVEVEL(matCENTER, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG, matCENTER);
q_ind = q_inds + matUINF(1,:);
fcolor = sqrt(sum(q_ind.^2,2));
p = patch('Faces',matDVE(:,:,1),'Vertices',matVLST,'FaceVertexCData',fcolor,'LineWidth',2);
p.FaceColor = 'flat';

% patch('Faces',matDVE(:,:,1),'Vertices',matVLST,'FaceVertexCData',matVLSTCP,'FaceColor','interp','LineWidth',2);
colorbar;

grid on
box on
axis equal

hold on
q_inds = fcnSDVEVEL(matCENTER, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG, matCENTER);
q_ind = q_inds + matUINF(1,:);
quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3), 'g')
hold off

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

% granularity = .5;
% x = -3:granularity:3;
% y = -3:granularity:3;
% z = -3:granularity:3;
% [X,Y,Z] = meshgrid(x,y,z);
% 
% s_ind = fcnSDVEVEL([X(:) Y(:) Z(:)],valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG, matCENTER);
% q_ind = s_ind + repmat(vecUINF, length(s_ind(:,1)),1);
% Xq = reshape(q_ind(:,1), size(X));
% Yq = reshape(q_ind(:,2), size(Y));
% Zq = reshape(q_ind(:,3), size(Z));
% 
% [Xs,Ys,Zs] = meshgrid(-3,y,z);
% hold on
% streamline(X,Y,Z,Xq,Yq,Zq,Xs,Ys,Zs);
% % quiver3(X(:),Y(:),Z(:),q_ind(:,1),q_ind(:,2),q_ind(:,3));
% hold off

% hFig213 = figure(213);
% clf(213);
% 
% idx_cp = matVLST(:,2) == 0;
% 
% theta = acosd(-matVLST(idx_cp,1));
% 
% scatter(theta, matVLSTCP(idx_cp), 'ok');
% hold on
% load('sphere_velocities.mat')
% plot(sphere_velocities(:,1), sphere_velocities(:,2), '--sb');
% grid on
% box on
% axis tight
% 
% xlabel('Chordwise Location (m)', 'FontSize',15);
% ylabel('Pressure Coefficient','FontSize',15);

%%
% v_row = 1:size(matVLST,1);
% 
% ii = 1;
% iii = 1;
% for jj = 1:length(v_row)
%     
%     [~, b] = ismembertol(v_row(jj), matDVE,'OutputAllIndices',true);
% 
%     [elem(:,1) elem(:,2)] = ind2sub(size(matDVE), cell2mat(b));
% 
%     for i = 1:size(elem,1)
%         point = matPLEX(elem(i,2),:,elem(i,1));
%         vort = [2.*matCOEFF(elem(i,1),3).*point(:,1) + matCOEFF(elem(i,1),4) + matCOEFF(elem(i,1),5).*point(:,2), 2.*matCOEFF(elem(i,1),1).*point(:,2) + matCOEFF(elem(i,1),2) + matCOEFF(elem(i,1),5).*point(:,1), point(:,2).*0];
%         vort_glob(ii,:) = fcnSTARGLOB(vort, matROTANG(elem(i,1),1), matROTANG(elem(i,1),2), matROTANG(elem(i,1),3));
%         point_glob(ii,:) = matVLST(v_row(jj),:);
%         v_num(ii,:) = v_row(jj);
%         ii = ii + 1;
%     end
%     
%     vort_avg(iii,:) = mean(vort_glob(ii-size(elem,1):ii-1,:),1);
%     point_avg(iii,:) = point_glob(ii-1,:);
%     
%     iii = iii + 1;
%     elem = [];
%     point = [];
%     vort = [];
% end
% 
% hFig28 = figure(28);
% clf(28);
% 
% subplot(1,2,1)
% scatter(v_num, vort_glob(:,1));
% ylabel('Gamma_x','FontSize',15);
% xlabel('Vertex Number','FontSize',15);
% 
% grid minor
% box on
% axis tight
% 
% subplot(1,2,2)
% scatter(v_num, vort_glob(:,2));
% ylabel('Gamma_y','FontSize',15);
% xlabel('Vertex Number','FontSize',15);
% 
% grid minor
% box on
% axis tight

