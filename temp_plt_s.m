q1 = q(:,:,1);
q2 = q(:,:,2);
q3 = q(:,:,3);

% GLOBAL ----------------------------------------------------------------------------------------------------------------------------
hFig1 = figure(28);
clf(28);

ax(1) = subplot(1,2,1);
patch('Faces',matDVE(:,:,1),'Vertices',matVLST,'FaceColor','r','LineWidth',2);
alpha(0);
hold on

for ii = 1:length(matVLST(:,1))
    str = sprintf('P%d',ii);
    text(matVLST(ii,1),matVLST(ii,2),matVLST(ii,3),str,'Color','r','FontSize',20);
end

scatter3(fpg(:,1), fpg(:,2), fpg(:,3), 100, 'ok', 'filled')
pstr = sprintf('PA (%0.2g,%0.2g,%0.2g)',fpg(1),fpg(2),fpg(3));
text(fpg(:,1), fpg(:,2), fpg(:,3), pstr, 'FontSize',15,'Color','m')

hold off
axis equal
box off
grid minor

title('Global Coordinate System','FontSize',15);
xlabel('X', 'FontSize', 15);
ylabel('Y', 'FontSize', 15);
zlabel('Z', 'FontSize', 15);

% LOCAL ------------------------------------------------------------------------------------------------------------------------------
ax(2) = subplot(1,2,2);
patch(matPLEX(:,1), matPLEX(:,2), matPLEX(:,3),'LineWidth',2,'EdgeColor','k','FaceAlpha',0);

hold on
scatter3(pa(:,1), pa(:,2), pa(:,3), 100, 'ok', 'filled')
pstr = sprintf('PA (%0.2g,%0.2g,%0.2g)',pa(1),pa(2),pa(3));
text(pa(:,1), pa(:,2), pa(:,3), pstr, 'FontSize',15,'Color','m')

% h1 = quiver3(0, 0, 0, matDVECT(:,1,1), matDVECT(:,2,1), matDVECT(:,3,1), 0.25, 'k'); % xsi
% h2 = quiver3(0, 0, 0, matDVECT(:,1,2), matDVECT(:,2,2), matDVECT(:,3,2), 0.25, 'b'); % eta
% h3 = quiver3(0, 0, 0, matDVECT(:,1,3), matDVECT(:,2,3), matDVECT(:,3,3), 0.25, 'm'); % zeta (normal)

%     scatter3(pb(:,1), pb(:,2), pb(:,3), 100, 'og', 'filled')
%     text(pb(:,1), pb(:,2), pb(:,3), 'PB', 'FontSize',20,'Color','m')
%
%     % (b1,b2,b3)
%     quiver3(pb(:,1),pb(:,2),pb(:,3), b1(1,1), b1(1,2), b1(1,3), 0.2, 'b', 'LineWidth',2)
%     quiver3(pb(:,1),pb(:,2),pb(:,3), b2(1,1), b2(1,2), b2(1,3), 0.2, 'b', 'LineWidth',2)
%     quiver3(pb(:,1),pb(:,2),pb(:,3), b3(1,1), b3(1,2), b3(1,3), 0.2, 'b', 'LineWidth',2)

for i = 1:3
    text(matPLEX(i,1), matPLEX(i,2), matPLEX(i,3), ['P' num2str(i)], 'FontSize',20,'Color','r')
end

hold off
axis equal
box off
grid minor

title('Local Coordinate System','FontSize',15);
xlabel('\xi', 'FontSize', 15);
ylabel('\eta', 'FontSize', 15);
zlabel('\zeta', 'FontSize', 15);

view([0 90])

%     linkprop(ax,{'xlim','ylim','zlim','CameraUpVector', 'CameraPosition', 'CameraTarget'});

% LOCAL AGAIN ----------------------------------------------------------------------------------
fsiz = 15;
hFig57 = figure(57);
clf(57);

ax2(1) = subplot(2,2,1);
patch(matPLEX(:,1), matPLEX(:,2), matPLEX(:,3),'LineWidth',2,'EdgeColor','k','FaceAlpha',0,'LineStyle','--');

hold on
%     scatter3(pa(:,1), pa(:,2), pa(:,3), 100, 'ok', 'filled')
%     text(pa(:,1), pa(:,2), pa(:,3), 'PA', 'FontSize',fsiz,'Color','m')

scatter3(pb(:,1), pb(:,2), pb(:,3), 100, 'og', 'filled')
% text(pb(:,1), pb(:,2), pb(:,3), 'PB', 'FontSize',fsiz,'Color','m')

% (b1,b2,b3)
quiver3(pb(:,1),pb(:,2),pb(:,3), b(1,1,1), b(1,2,1), b(1,3,1), 0.2, 'b', 'LineWidth',1)
text(pb(:,1)+0.2*b(1,1,1), pb(:,2)+0.2*b(1,2,1), pb(:,3)+0.2*b(1,3,1), 'B1','FontSize',fsiz,'Color','b')
quiver3(pb(:,1),pb(:,2),pb(:,3), b(1,1,2), b(1,2,2), b(1,3,2), 0.2, 'b', 'LineWidth',1)
text(pb(:,1)+0.2*b(1,1,2), pb(:,2)+0.2*b(1,2,2), pb(:,3)+0.2*b(1,3,2), 'B2','FontSize',fsiz,'Color','b')
quiver3(pb(:,1),pb(:,2),pb(:,3), b(1,1,3), b(1,2,3), b(1,3,3), 0.2, 'b', 'LineWidth',1)
text(pb(:,1)+0.2*b(1,1,3), pb(:,2)+0.2*b(1,2,3), pb(:,3)+0.2*b(1,3,3), 'B3','FontSize',fsiz,'Color','b')

for i = 1:3
    text(matPLEX(i,1), matPLEX(i,2), matPLEX(i,3), ['P' num2str(i)], 'FontSize',fsiz,'Color','r')
end

hold off
axis equal
box off
grid minor

xlabel('\xi', 'FontSize', 15);
ylabel('\eta', 'FontSize', 15);
zlabel('\zeta', 'FontSize', 15);

view([0 90])

clear q1 q2 q3