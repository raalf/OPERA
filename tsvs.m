clc
clear

endpoints(:,:,1) = [0 0 0];
endpoints(:,:,2) = [1 1 0];

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

y = 0:granularity:1.5;
% z = 0

[X,Y,Z] = meshgrid(x,y,z);

fpl = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');

% fpl = [-4.5 -5 0];

test = 1;

len = length(fpl(:,1));

%% One sheet using old way
temp_endpoints = repmat(endpoints,len,1,1);
temp_phi = repmat(phi,len,1);
temp_yaw = repmat(yaw,len,1);
temp_k = repmat(k,len,1);

limits = [ones(len,1) -ones(len,1)];

hspan1 = abs(endpoints(:,1,2) - endpoints(:,1,1))./2;
hspan2 = abs(endpoints(:,2,2) - endpoints(:,2,1))./2;
hspan = hspan1.*cos(yaw) + hspan2.*sin(yaw);

hspan = repmat(hspan, len, 1);

[aloc, bloc, cloc] = fcnVSIND(temp_endpoints, temp_phi, temp_yaw, fpl, temp_k);

% [aloc, bloc, cloc] = fcnHORSTMANN(temp_endpoints, hspan, temp_phi, temp_yaw, fpl, temp_k);


%%


D = [cloc bloc aloc];
D = reshape(reshape(D', 1, 9, []), 3, 3, len);

q_ind = permute(sum(D.*repmat(reshape(coeff',1,3,[]),3,1,len),2),[2 1 3]);
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

