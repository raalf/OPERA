function [hFig1] = fcnPLOTBODY(verbose, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, vecUINF)
% This function plots all elements, and can label vertices, faces and edges.

hFig1 = figure(1);
clf(1);

patch('Faces',matDVE(:,:,1),'Vertices',matVLST,'FaceColor','r','LineWidth',2);
alpha(0);
hold on

if verbose == 1
    for ii = 1:valNELE
        str = sprintf('%d',ii);
        text(matCENTER(ii,1),matCENTER(ii,2),matCENTER(ii,3),str,'Color','k','FontSize',20);
    end
    
    for ii = 1:length(matVLST(:,1))
        str = sprintf('%d',ii);
        text(matVLST(ii,1),matVLST(ii,2),matVLST(ii,3),str,'Color','g','FontSize',20);
    end
    
    edge1 = matVLST(matELST(:,1),:);
    edge2 = matVLST(matELST(:,2),:);
    mid = (edge1+edge2)./2;
    for ii = 1:length(mid)
        str = sprintf('%d',ii);
        text(mid(ii,1),mid(ii,2),mid(ii,3),str,'Color','b','FontSize',20);
    end
    
    [q_ind] = fcnINDVEL(1:valNELE, matCENTER, matCOEFF, matDVE, matDVECT, matVLST, matPLEX);
    
    quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), q_ind(:,1)+vecUINF(1), q_ind(:,2)+vecUINF(2), q_ind(:,3)+vecUINF(3), 0.25, 'g')
    
    quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), matDVECT(:,1,1), matDVECT(:,2,1), matDVECT(:,3,1), 0.25, 'b') % eta
    quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), matDVECT(:,1,2), matDVECT(:,2,2), matDVECT(:,3,2), 0.25, 'k') % xi
    quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), matDVECT(:,1,3), matDVECT(:,2,3), matDVECT(:,3,3), 0.25,'m') % zeta (normal)
    
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

