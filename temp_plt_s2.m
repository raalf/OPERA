ax2(ii+1) = subplot(2,2,ii+1);
copyobj(get(ax2(1),'Children'),ax2(ii+1));

qm_l = fcnSTARGLOB(qm, loc_rot(:,1), loc_rot(:,2), loc_rot(:,3));
qn_l = fcnSTARGLOB(qn, loc_rot(:,1), loc_rot(:,2), loc_rot(:,3));
c1_l = fcnSTARGLOB(c1, loc_rot(:,1), loc_rot(:,2), loc_rot(:,3));
c2_l = fcnSTARGLOB(c2, loc_rot(:,1), loc_rot(:,2), loc_rot(:,3));
c3_l = fcnSTARGLOB(c3, loc_rot(:,1), loc_rot(:,2), loc_rot(:,3));

p_l = [pb; pb+qm_l; pb+qn_l];
patch(p_l(:,1), p_l(:,2), p_l(:,3),'LineWidth',2,'EdgeColor','k','FaceAlpha',0);

for i = 1:3
    text(p_l(i,1), p_l(i,2), p_l(i,3), ['P' num2str(i)], 'FontSize',20,'Color','r')
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
quiver3(pb(:,1), pb(:,2), pb(:,3), qm_l(:,1), qm_l(:,2), qm_l(:,3), 0, 'm','MaxHeadSize',mhs)
text(pb(:,1) + laa*qm_l(:,1), pb(:,2) + laa*qm_l(:,2), pb(:,3) + laa*qm_l(:,3), 'qm_l', 'FontSize',fsiz)

quiver3(pb(:,1), pb(:,2), pb(:,3), qn_l(:,1), qn_l(:,2), qn_l(:,3), 0, 'm','MaxHeadSize',mhs)
text(pb(:,1) + laa*qn_l(:,1), pb(:,2) + laa*qn_l(:,2), pb(:,3) + laa*qn_l(:,3), 'qn_l', 'FontSize',fsiz)

quiver3(pb(:,1),pb(:,2),pb(:,3), c1_l(1,1), c1_l(1,2), c1_l(1,3), 0.2, 'r', 'LineWidth',1)
text(pb(:,1)+0.3*c1_l(1,1), pb(:,2)+0.3*c1_l(1,2), pb(:,3)+0.3*c1_l(1,3), 'C1','FontSize',fsiz,'Color','r')
quiver3(pb(:,1),pb(:,2),pb(:,3), c2_l(1,1), c2_l(1,2), c2_l(1,3), 0.2, 'r', 'LineWidth',1)
text(pb(:,1)+0.3*c2_l(1,1), pb(:,2)+0.3*c2_l(1,2), pb(:,3)+0.3*c2_l(1,3), 'C2','FontSize',fsiz,'Color','r')
quiver3(pb(:,1),pb(:,2),pb(:,3), c3_l(1,1), c3_l(1,2), c3_l(1,3), 0.2, 'r', 'LineWidth',1)
% text(pb(:,1)+0.3*c3_l(1,1), pb(:,2)+0.3*c3_l(1,2), pb(:,3)+0.3*c3_l(1,3), 'C3','FontSize',fsiz,'Color','r')


hold off

strtitle = sprintf('S%d (pb, p%d, p%d)', ii, ii, ii+1);
title(strtitle,'FontSize',15)
xlabel('\xi', 'FontSize', 15);
ylabel('\eta', 'FontSize', 15);
zlabel('\zeta', 'FontSize', 15);