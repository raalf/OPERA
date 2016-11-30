clc
clear

endpoints(:,:,1) = [-1 0 0];
endpoints(:,:,2) = [1 0 0];

phi = 0;
yaw = 0;

k = 0.1;

granularity = 0.5;

x = -2:granularity:2;
y = -2:granularity:2;
z = -2:granularity:2;

[X,Y,Z] = meshgrid(x,y,z);

fpl = [reshape(X, [], 1, 1) reshape(Y, [], 1, 1) reshape(Z, [], 1, 1)];

len = length(fpl(:,1));

temp_endpoints = repmat(endpoints,len,1,1);
temp_phi = repmat(phi,len,1);
temp_yaw = repmat(yaw,len,1);
temp_k = repmat(k,len,1);

lim1 = -1;
lim2 = 1;

[aloc, bloc, cloc] = fcnHVSIND(temp_endpoints, lim1, lim2, temp_phi, temp_yaw, fpl, temp_k);

D = [cloc bloc aloc];
D = reshape(reshape(D', 1, 9, []), 3, 3, len);

q_ind = permute(sum(D.*repmat(reshape([0 1 0]',1,3,[]),3,1,1),2),[2 1 3]);
q_ind = reshape(permute(q_ind,[3 1 2]),[],3,1)./(-4*pi);

% Plot
hFig12 = figure(12);
clf(12);

plot3(reshape(endpoints(:,1,:),[],1,1), reshape(endpoints(:,2,:),[],1,1), reshape(endpoints(:,3,:),[],1,1),'-k','LineWidth',2)
hold on
quiver3(fpl(:,1), fpl(:,2), fpl(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3));
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
clc
clear

% First sheet
endpoints(:,:,1) = [-1 0 0];
endpoints(:,:,2) = [1 0 0];

phi = 0;
yaw = 0;

k = 0.1;

granularity = 0.5;

x = -2:granularity:2;
y = -2:granularity:2;
z = -2:granularity:2;

[X,Y,Z] = meshgrid(x,y,z);

fpl = [reshape(X, [], 1, 1) reshape(Y, [], 1, 1) reshape(Z, [], 1, 1)];

len = length(fpl(:,1));

temp_endpoints = repmat(endpoints,len,1,1);
temp_phi = repmat(phi,len,1);
temp_yaw = repmat(yaw,len,1);
temp_k = repmat(k,len,1);

lim1 = -1;
lim2 = 0;

[aloc, bloc, cloc] = fcnHVSIND(temp_endpoints, lim1, lim2, temp_phi, temp_yaw, fpl, temp_k);

D = [cloc bloc aloc];
D = reshape(reshape(D', 1, 9, []), 3, 3, len);

q_ind1 = permute(sum(D.*repmat(reshape([0 1 0]',1,3,[]),3,1,1),2),[2 1 3]);
q_ind1 = reshape(permute(q_ind1,[3 1 2]),[],3,1)./(-4*pi);

% Second sheet
endpoints(:,:,1) = [-1 0 0];
endpoints(:,:,2) = [1 0 0];

phi = 0;
yaw = 0;

k = 0.1;

granularity = 0.5;

x = -2:granularity:2;
y = -2:granularity:2;
z = -2:granularity:2;

[X,Y,Z] = meshgrid(x,y,z);

fpl = [reshape(X, [], 1, 1) reshape(Y, [], 1, 1) reshape(Z, [], 1, 1)];

len = length(fpl(:,1));

temp_endpoints = repmat(endpoints,len,1,1);
temp_phi = repmat(phi,len,1);
temp_yaw = repmat(yaw,len,1);
temp_k = repmat(k,len,1);

lim1 = 0;
lim2 = 1;

[aloc, bloc, cloc] = fcnHVSIND(temp_endpoints, lim1, lim2, temp_phi, temp_yaw, fpl, temp_k);

D = [cloc bloc aloc];
D = reshape(reshape(D', 1, 9, []), 3, 3, len);

q_ind2 = permute(sum(D.*repmat(reshape([0 1 0]',1,3,[]),3,1,1),2),[2 1 3]);
q_ind2 = reshape(permute(q_ind2,[3 1 2]),[],3,1)./(-4*pi);

q_ind = q_ind1 + q_ind2;

% Plot
hFig13 = figure(13);
clf(13);

plot3(reshape(endpoints(:,1,:),[],1,1), reshape(endpoints(:,2,:),[],1,1), reshape(endpoints(:,3,:),[],1,1),'-k','LineWidth',2)
hold on
quiver3(fpl(:,1), fpl(:,2), fpl(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3));
hold off

view([0 0]);

grid on
axis tight
axis equal
box on


xlabel('X-Dir','FontSize',15);
ylabel('Y-Dir','FontSize',15);
zlabel('Z-Dir','FontSize',15);