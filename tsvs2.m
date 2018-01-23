clc
clear


granularity = .05;
% x = -0.5:granularity:1.2;
y = -0.2:granularity:1.2;
% y = ones(size(x)) - 0.5;
z = -0.2:granularity:0.2;

% z = 0.01;
x = 0:granularity:1;
% x = zeros(size(z)) + 0.33;
% y = zeros(size(z)) + 0.5;


[X,Y,Z] = meshgrid(x,y,z);
fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');
fpl = fpg;



len = size(fpl,1);
coeff = [0 1 0];

temp_yaw = repmat(pi/2,len,1);
temp_k = repmat(0.01,len,1);


hFig13 = figure(13);
clf(13);

%% Sheet 1
endpoints(:,:,1) = [0 0 0];
endpoints(:,:,2) = [1 0.5 0];
phi = -atan((endpoints(:,1,2)-endpoints(:,1,1))./(endpoints(:,2,2)-endpoints(:,2,1)));

temp_endpoints = repmat(endpoints, len, 1);
phi = repmat(phi, len, 1);

[aloc, bloc, cloc] = fcnVSIND(temp_endpoints, phi, temp_yaw, fpl, temp_k);

D = [cloc bloc aloc];
D = reshape(reshape(D', 1, 9, []), 3, 3, len);

q_ind1 = permute(sum(D.*repmat(reshape(coeff',1,3,[]),3,1,1),2),[2 1 3]);
q_ind1 = reshape(permute(q_ind1,[3 1 2]),[],3,1)./(-4*pi);

plot3(reshape(endpoints(:,1,:),[],1,1), reshape(endpoints(:,2,:),[],1,1), reshape(endpoints(:,3,:),[],1,1),'-k','LineWidth',2)


%% Second sheet
endpoints = [];
endpoints(:,:,2) = [0 1 0];
endpoints(:,:,1) = [1 0.5 0];
phi = -atan((endpoints(:,1,2)-endpoints(:,1,1))./(endpoints(:,2,2)-endpoints(:,2,1)));

temp_endpoints = repmat(endpoints, len, 1);
phi = repmat(phi, len, 1);

[aloc, bloc, cloc] = fcnVSIND(temp_endpoints, phi, temp_yaw, fpl, temp_k);

D = [cloc bloc aloc];
D = reshape(reshape(D', 1, 9, []), 3, 3, len);

q_ind2 = permute(sum(D.*repmat(reshape(coeff',1,3,[]),3,1,1),2),[2 1 3]);
q_ind2 = reshape(permute(q_ind2,[3 1 2]),[],3,1)./(-4*pi);

% q_ind = q_ind2;
q_ind = q_ind1 - q_ind2;

hold on
plot3(reshape(endpoints(:,1,:),[],1,1), reshape(endpoints(:,2,:),[],1,1), reshape(endpoints(:,3,:),[],1,1),'-k','LineWidth',2)


%% Plot

quiver3(fpl(:,1), fpl(:,2), fpl(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3),'b');

hold off


grid on
axis tight
axis equal
box on


xlabel('X-Dir','FontSize',15);
ylabel('Y-Dir','FontSize',15);
zlabel('Z-Dir','FontSize',15);
