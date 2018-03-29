clc
clear

load('matlab.mat')

len = size(matUINF,1);
opera_sphere_vel = sqrt(sum((q_inds + matUINF).^2,2))./sqrt(sum(matUINF.^2,2));
vec_prj = (matCENTER - (dot(matCENTER, repmat([0 1 0], len, 1),2)).*repmat([0 1 0], len, 1));
vec_prj = vec_prj./(sqrt(sum(vec_prj.^2,2)));
theta = acosd(dot(vec_prj, repmat([1 0 0], len, 1),2));

hFig213 = figure(213);
clf(213);
hold on
load('panair_sphere_vel.mat');
plot(panair_sphere_vel(:,1), panair_sphere_vel(:,2), '--b');
scatter(theta, opera_sphere_vel, 'ok');
grid minor
box on
axis tight
xlabel('Theta (Degrees)', 'FontSize',15);
ylabel('Velocity Ratio (V/V_{inf})','FontSize',15);