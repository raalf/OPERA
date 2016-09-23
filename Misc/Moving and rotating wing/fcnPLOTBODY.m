function [hFig1] = fcnPLOTBODY(verbose, DVE, NELE, VLST, ELST, DVECT, CENTER)
% This function plots all elements, and can label vertices, faces and edges.

hFig1 = figure(1);
% clf(1);

TR = triangulation(DVE(:,:,1),VLST);
trisurf(TR,'FaceColor',[.5 .5 .5],'LineWidth',2);
alpha(0);
hold on

if verbose == 1
    for ii = 1:NELE
        str = sprintf('%d',ii);
        text(CENTER(ii,1),CENTER(ii,2),CENTER(ii,3),str,'Color','k','FontSize',20);
    end
    
    for ii = 1:length(VLST(:,1))
        str = sprintf('%d',ii);
        text(VLST(ii,1),VLST(ii,2),VLST(ii,3),str,'Color','g','FontSize',20);
    end
    
    edge1 = VLST(ELST(:,1),:);
    edge2 = VLST(ELST(:,2),:);
    mid = (edge1+edge2)./2;
    for ii = 1:length(mid)
        str = sprintf('%d',ii);
        text(mid(ii,1),mid(ii,2),mid(ii,3),str,'Color','b','FontSize',20);
    end
    
    quiver3(CENTER(:,1),CENTER(:,2),CENTER(:,3), DVECT(:,1,1), DVECT(:,2,1), DVECT(:,3,1), 0.25, 'b') % eta
    quiver3(CENTER(:,1),CENTER(:,2),CENTER(:,3), DVECT(:,1,2), DVECT(:,2,2), DVECT(:,3,2), 0.25, 'k') % xi
    quiver3(CENTER(:,1),CENTER(:,2),CENTER(:,3), DVECT(:,1,3), DVECT(:,2,3), DVECT(:,3,3), 0.25,'m') % zeta (normal)
    
end

hold off
set(gcf,'Renderer','OpenGL');
axis equal
axis tight
box on
grid on

xlabel('X-Dir', 'FontSize', 15);
ylabel('Y-Dir', 'FontSize', 15);
zlabel('Z-Dir', 'FontSize', 15);

end

