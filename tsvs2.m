clc
clear

endpoints(1,:,1) = [0 0 0];
endpoints(1,:,2) = [1 1 0];

endpoints(2,:,1) = [0 0 0];
endpoints(2,:,2) = [1 0 0];

coeff = [1 1 0];

phi = atan((endpoints(:,2,2)-endpoints(:,2,1))./(endpoints(:,1,2)-endpoints(:,1,1)));
yaw = 0;

k = 0.1;

granularity = 0.25;
xlims = 5;
ylims = 5;
zlims = 2;

x = -xlims:granularity:xlims;
y = -ylims:granularity:ylims;
z = -zlims:granularity:zlims;

y = -.5:granularity:1.5;
% z = 0

[X,Y,Z] = meshgrid(x,y,z);

fpl = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');

% fpl = [0.5 0.5 0];
% fpl = [1.25 1.25 0]
fpl = [0.75 0.5 0]

test = 1;

len = length(fpl(:,1));

% One sheet using old way
temp_endpoints = repmat(endpoints(1,:,:),len,1,1);
temp_phi = repmat(phi(1),len,1);
temp_yaw = repmat(yaw,len,1);
temp_k = repmat(k,len,1);

hspan1 = abs(endpoints(1,1,2) - endpoints(1,1,1))./2;
hspan2 = abs(endpoints(1,2,2) - endpoints(1,2,1))./2;
hspan = hspan1.*cos(yaw) + hspan2.*sin(yaw);

hspan = repmat(hspan, len, 1);

[aloc1, bloc1, cloc1] = fcnVSIND(temp_endpoints, temp_phi, temp_yaw, fpl, temp_k);

% [aloc, bloc, cloc] = fcnHORSTMANN(temp_endpoints, hspan, temp_phi, temp_yaw, fpl, temp_k);

% Second sheet using old way
temp_endpoints = repmat(endpoints(2,:,:),len,1,1);
temp_phi = repmat(phi(2),len,1);

hspan1 = abs(endpoints(2,1,2) - endpoints(2,1,1))./2;
hspan2 = abs(endpoints(2,2,2) - endpoints(2,2,1))./2;
hspan = hspan1.*cos(yaw) + hspan2.*sin(yaw);

hspan = repmat(hspan, len, 1);

[aloc2, bloc2, cloc2] = fcnVSIND(temp_endpoints, temp_phi, temp_yaw, fpl, temp_k);

%%
% u = [0.866 0 .5]

aloc = aloc2 - aloc1;
bloc = bloc2 - bloc1;
cloc = cloc2 - cloc1;


D = [cloc bloc aloc];
D = reshape(reshape(D', 1, 9, []), 3, 3, len);

q_ind = permute(sum(D.*repmat(reshape(coeff',1,3,[]),3,1,len),2),[2 1 3]);
q_ind = reshape(permute(q_ind,[3 1 2]),[],3,1)./(-4*pi);

q_ind = q_ind;

% Plot
hFig12 = figure(12);
clf(12);

plot3(reshape(endpoints(1,1,:),[],1,1), reshape(endpoints(1,2,:),[],1,1), reshape(endpoints(1,3,:),[],1,1),'-k','LineWidth',2)
hold on
plot3(reshape(endpoints(2,1,:),[],1,1), reshape(endpoints(2,2,:),[],1,1), reshape(endpoints(2,3,:),[],1,1),'-k','LineWidth',2)
hold off
hold on
quiver3(fpl(:,1), fpl(:,2), fpl(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3),'b');
hold off

view([0 0]);

grid on
axis tight
axis equal
box on

xlabel('X-Dir','FontSize',15);
ylabel('Y-Dir','FontSize',15);
zlabel('Z-Dir','FontSize',15);

%%
hold on
D = [cloc1 bloc1 aloc1];
D = reshape(reshape(D', 1, 9, []), 3, 3, len);

q_ind = permute(sum(D.*repmat(reshape(coeff',1,3,[]),3,1,len),2),[2 1 3]);
q_ind = reshape(permute(q_ind,[3 1 2]),[],3,1)./(-4*pi);

quiver3(fpl(:,1), fpl(:,2), fpl(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3),0,'g');

D = [cloc2 bloc2 aloc2];
D = reshape(reshape(D', 1, 9, []), 3, 3, len);

q_ind = permute(sum(D.*repmat(reshape(coeff',1,3,[]),3,1,len),2),[2 1 3]);
q_ind = reshape(permute(q_ind,[3 1 2]),[],3,1)./(-4*pi);

quiver3(fpl(:,1), fpl(:,2), fpl(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3),0,'m');
hold off


