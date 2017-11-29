clc
clear

p0 = [2 -1; 1 0; 5 2];

% Edge vectors by going around in order
edge = [p0(1,:) - p0(2,:); p0(2,:) - p0(3,:); p0(3,:) - p0(1,:)];

% Finding normals and normalizing
edge_normal = [edge(:,2) -edge(:,1)];
edge_normal = edge_normal./(sqrt(sum(edge_normal(:,1).^2 + edge_normal(:,2).^2,2)));

edge_midpoints = [mean(p0(1:2,:),1); mean(p0(2:3,:),1); mean(p0(3:-2:1,:),1)];

check_normal = dot(edge_normal(1,:), edge(3,:));
check_normal(check_normal > 0) = 1;

if check_normal == 1
   edge_normal = [-edge_normal(:,1) -edge_normal(:,2)]; 
end

%%
hfig15 = figure(15);
clf(15)
fill(p0(:,1), p0(:,2), 'w', 'FaceAlpha',0)
grid minor
box on
axis equal

hold on
quiver(edge_midpoints(:,1), edge_midpoints(:,2), edge_normal(:,1), edge_normal(:,2))
hold off