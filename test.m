clc
clear

load('temp1.mat')

%%
q_inds = fcnSDVEVEL(matVLST, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG, matCENTER);
q_ind = q_inds + matUINF(1,:);

matVLSTCP = sqrt(sum(q_ind.^2,2));

hFig1 = figure(1);
clf(1);

patch('Faces',matDVE(:,:,1),'Vertices',matVLST,'FaceVertexCData',matVLSTCP,'FaceColor','interp','LineWidth',2);
colorbar;

grid on
box on
axis equal

% hold on

% quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3), 'g')
% hold off

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

fpg = [-0.25 0 3];

[s_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG, matCENTER);

% q_ind = s_ind + repmat(vecUINF, length(s_ind(:,1)),1);
q_ind = s_ind;
hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3))
hold off
axis tight