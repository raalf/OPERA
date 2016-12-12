clc
clear

endpoints(:,:,1) = [-0.5 0 0];
endpoints(:,:,2) = [0.5 0 0];

phi = atan((endpoints(:,2,2)-endpoints(:,2,1))./(endpoints(:,1,2)-endpoints(:,1,1)));
yaw = 0;

k = 0.1;

granularity = 0.5;
xlims = 10;
ylims = 10;
zlims = 10;

x = -xlims:granularity:xlims;
y = -ylims:granularity:ylims;
z = -zlims:granularity:zlims;

[X,Y,Z] = meshgrid(x,y,z);

fpl = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');

% fpl = [-5 0 1];

test = 1;

len = length(fpl(:,1));

%% One sheet using old way
temp_endpoints = repmat(endpoints,len,1,1);
temp_phi = repmat(phi,len,1);
temp_yaw = repmat(yaw,len,1);
temp_k = repmat(k,len,1);

limits = [ones(len*5,1) -ones(len*5,1)];

hspan1 = abs(endpoints(:,1,2) - endpoints(:,1,1))./2;
hspan2 = abs(endpoints(:,2,2) - endpoints(:,2,1))./2;
hspan = hspan1.*cos(yaw) + hspan2.*sin(yaw);

hspan = repmat(hspan, len, 1);

[aloc, bloc, cloc] = fcnHVSIND2(temp_endpoints, [-1 1], hspan, temp_phi, temp_yaw, fpl, temp_k);

D = [cloc bloc aloc];
D = reshape(reshape(D', 1, 9, []), 3, 3, len);

q_ind = permute(sum(D.*repmat(reshape([0 1 0]',1,3,[]),3,1,len),2),[2 1 3]);
q_ind = reshape(permute(q_ind,[3 1 2]),[],3,1)./(-4*pi);

q_ind(test,:)

% Plot
hFig12 = figure(12);
clf(12);

plot3(reshape(endpoints(:,1,:),[],1,1), reshape(endpoints(:,2,:),[],1,1), reshape(endpoints(:,3,:),[],1,1),'-k','LineWidth',2)
hold on
quiver3(fpl(:,1), fpl(:,2), fpl(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3),0,'b');
hold off

view([0 0]);

grid on
axis tight
axis equal
box on

xlabel('X-Dir','FontSize',15);
ylabel('Y-Dir','FontSize',15);
zlabel('Z-Dir','FontSize',15);

%% Adding two sheets, the new way

% midpoint = [0.5 0.5 0];
% 
% % First Sheet
% limits = repmat([1 0],len,1);
% 
% % Half-spans of all elements
% hspan1 = abs(temp_endpoints(:,1,2) - temp_endpoints(:,1,1))./2;
% hspan2 = abs(temp_endpoints(:,2,2) - temp_endpoints(:,2,1))./2;
% hspan = hspan1.*cos(temp_yaw) + hspan2.*sin(temp_yaw);
% 
% [aloc, bloc, cloc] = fcnHVSIND(temp_endpoints, limits, hspan, temp_phi, temp_yaw, fpl, temp_k);
% 
% D = [cloc bloc aloc];
% D = reshape(reshape(D', 1, 9, []), 3, 3, len);
% 
% q_ind1 = permute(sum(D.*repmat(reshape([0 1 0]',1,3,[]),3,1,1),2),[2 1 3]);
% q_ind1 = reshape(permute(q_ind1,[3 1 2]),[],3,1)./(-4*pi);
% 
% % Second sheet
% limits = repmat([0 -1],len,1);
% 
% % Half-spans of all elements
% hspan1 = abs(temp_endpoints(:,1,2) - temp_endpoints(:,1,1))./2;
% hspan2 = abs(temp_endpoints(:,2,2) - temp_endpoints(:,2,1))./2;
% hspan = hspan1.*cos(temp_yaw) + hspan2.*sin(temp_yaw);
% 
% [aloc, bloc, cloc] = fcnHVSIND(temp_endpoints, limits, hspan, temp_phi, temp_yaw, fpl, temp_k);
% 
% D = [cloc bloc aloc];
% D = reshape(reshape(D', 1, 9, []), 3, 3, len);
% 
% q_ind2 = permute(sum(D.*repmat(reshape([0 1 0]',1,3,[]),3,1,1),2),[2 1 3]);
% q_ind2 = reshape(permute(q_ind2,[3 1 2]),[],3,1)./(-4*pi);
% 
% q_ind = q_ind1 + q_ind2;
% 
% q_ind(test,:)
% 
% % Plot
% hFig13 = figure(13);
% clf(13);
% 
% plot3(reshape(endpoints(:,1,:),[],1,1), reshape(endpoints(:,2,:),[],1,1), reshape(endpoints(:,3,:),[],1,1),'-k','LineWidth',2)
% hold on
% scatter3(midpoint(:,1), midpoint(:,2), midpoint(:,3),100,'filled','r');
% quiver3(fpl(:,1), fpl(:,2), fpl(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3),'b');
% hold off
% 
% view([0 0]);
% 
% grid on
% axis tight
% axis equal
% box on
% 
% 
% xlabel('X-Dir','FontSize',15);
% ylabel('Y-Dir','FontSize',15);
% zlabel('Z-Dir','FontSize',15);