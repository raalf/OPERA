clc
clear

load('opera_sphere_25.mat')
%% Plot with Streamlines
hFig1 = figure(1);
clf(1);
hold on
q_ind = q_inds + matUINF(1,:);
fcolor = sqrt(sum(q_ind.^2,2));
p = patch('Faces',matDVE(:,:,1),'Vertices',matVLST,'FaceVertexCData',fcolor,'LineWidth',0.1);
p.FaceColor = 'flat';
% quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3), 'g')
hold off
colorbar;
grid on
box on
axis equal

granularity = .5;
x = -1.5:granularity:1.5;
% y = -3:granularity:3;
y = -1.5:granularity:1.5;
z = -1.5:granularity:1.5;
[X,Y,Z] = meshgrid(x,y,z);

q_ind = s_ind + repmat(matUINF(1,:), length(s_ind(:,1)),1);
Xq = reshape(q_ind(:,1), size(X));
Yq = reshape(q_ind(:,2), size(Y));
Zq = reshape(q_ind(:,3), size(Z));

[Xs,Ys,Zs] = meshgrid(-1.5,y,z);
hold on
streamline(X,Y,Z,Xq,Yq,Zq,Xs,Ys,Zs);
% quiver3(X(:),Y(:),Z(:),q_ind(:,1),q_ind(:,2),q_ind(:,3));
hold off

axis tight

%%
% len = size(matUINF,1);
% opera_sphere_vel = sqrt(sum((q_inds + matUINF).^2,2))./sqrt(sum(matUINF.^2,2));
% vec_prj = (matCENTER - (dot(matCENTER, repmat([0 1 0], len, 1),2)).*repmat([0 1 0], len, 1));
% vec_prj = vec_prj./(sqrt(sum(vec_prj.^2,2)));
% theta = acosd(dot(vec_prj, repmat([1 0 0], len, 1),2));
% 
% hFig213 = figure(213);
% clf(213);
% hold on
% load('panair_sphere_vel.mat');
% plot(panair_sphere_vel(:,1), panair_sphere_vel(:,2), '--b');
% scatter(theta, opera_sphere_vel, 'ok');
% grid minor
% box on
% axis tight
% xlabel('Theta (Degrees)', 'FontSize',15);
% ylabel('Velocity Ratio (V/V_{inf})','FontSize',15);

