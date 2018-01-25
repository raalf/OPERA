ax2(ii+1) = subplot(2,2,ii+1);
% copyobj(get(ax2(1),'Children'),ax2(ii+1));

p = [pb; pb+qm; pb+qn];
patch(p(:,1), p(:,2), p(:,3),'LineWidth',2,'EdgeColor','k','FaceAlpha',0);

for i = 1:3
    text(p(i,1), p(i,2), p(i,3), ['P' num2str(i)], 'FontSize',20,'Color','r')
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

fsiz = 15;
mhs = 0.5;
laa = 0.8;

hold on
quiver3(pb(:,1), pb(:,2), pb(:,3), qm(:,1), qm(:,2), qm(:,3), 0, 'm','MaxHeadSize',mhs)
text(pb(:,1) + laa*qm(:,1), pb(:,2) + laa*qm(:,2), pb(:,3) + laa*qm(:,3), 'qm', 'FontSize',fsiz)

quiver3(pb(:,1), pb(:,2), pb(:,3), qn(:,1), qn(:,2), qn(:,3), 0, 'm','MaxHeadSize',mhs)
text(pb(:,1) + laa*qn(:,1), pb(:,2) + laa*qn(:,2), pb(:,3) + laa*qn(:,3), 'qn', 'FontSize',fsiz)

c2 = (qn - qm);
c2 = c2./sqrt(sum(c2.^2,2));
c1 = cross(c2,c3);

quiver3(pb(:,1),pb(:,2),pb(:,3), c1(1,1), c1(1,2), c1(1,3), 0.2, 'r', 'LineWidth',1)
text(pb(:,1)+0.3*c1(1,1), pb(:,2)+0.3*c1(1,2), pb(:,3)+0.3*c1(1,3), 'C1','FontSize',fsiz,'Color','r')
quiver3(pb(:,1),pb(:,2),pb(:,3), c2(1,1), c2(1,2), c2(1,3), 0.2, 'r', 'LineWidth',1)
text(pb(:,1)+0.3*c2(1,1), pb(:,2)+0.3*c2(1,2), pb(:,3)+0.3*c2(1,3), 'C2','FontSize',fsiz,'Color','r')
quiver3(pb(:,1),pb(:,2),pb(:,3), c3(1,1), c3(1,2), c3(1,3), 0.2, 'r', 'LineWidth',1)
% text(pb(:,1)+0.3*c3(1,1), pb(:,2)+0.3*c3(1,2), pb(:,3)+0.3*c3(1,3), 'C3','FontSize',fsiz,'Color','r')


hold off
axis equal
box off
grid minor


strtitle = sprintf('S%d (pb, p%d, p%d)', ii, ii, ii+1);
title(strtitle,'FontSize',15)
xlabel('\xi', 'FontSize', 15);
ylabel('\eta', 'FontSize', 15);
zlabel('\zeta', 'FontSize', 15);