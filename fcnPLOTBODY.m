function [hFig1] = fcnPLOTBODY(verbose, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, fig_size, rendertype)
% This function plots all elements, and can label vertices, faces and edges.

hFig1 = figure(1);
clf(1);

patch('Faces',matDVE(:,:,1),'Vertices',matVLST,'FaceColor','r','LineWidth',2);
alpha(0);
hold on

% edge1 = matVLST(matELST(:,1),:);
% edge2 = matVLST(matELST(:,2),:);
% mid = (edge1+edge2)./2;
% for ii = 1:length(mid)
%     str = sprintf('%d',ii);
%     text(mid(ii,1),mid(ii,2),mid(ii,3),str,'Color','b','FontSize',15);
% end

if verbose == 1
    for ii = 1:valNELE
        str = sprintf('%d',ii);
        text(matCENTER(ii,1),matCENTER(ii,2),matCENTER(ii,3),str,'Color','k','FontSize',20);
    end
    
    for ii = 1:length(matVLST(:,1))
        str = sprintf('%d',ii);
        text(matVLST(ii,1),matVLST(ii,2),matVLST(ii,3),str,'Color','r','FontSize',20);
    end
    
    edge1 = matVLST(matELST(:,1),:);
    edge2 = matVLST(matELST(:,2),:);
    mid = (edge1+edge2)./2;
    for ii = 1:length(mid)
        str = sprintf('%d',ii);
        text(mid(ii,1),mid(ii,2),mid(ii,3),str,'Color','b','FontSize',20);
    end


    len = length(matCENTER(:,1));
    %     quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3), 0.25, 'g')
%     quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), repmat(vecUINF(1),len,1), repmat(vecUINF(2),len,1), repmat(vecUINF(3),len,1), 0.25, 'r')
    
    h1 = quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), matDVECT(:,1,1), matDVECT(:,2,1), matDVECT(:,3,1), 0.25, 'k'); % xsi
    h2 = quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), matDVECT(:,1,2), matDVECT(:,2,2), matDVECT(:,3,2), 0.25, 'b'); % eta
    h3 = quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), matDVECT(:,1,3), matDVECT(:,2,3), matDVECT(:,3,3), 0.25, 'm'); % zeta (normal)
%     legend([h1,h2,h3],'Local Eta direction','Local Xsi direction','Normal','Location','NorthWest')
    
end

hold off
set(gcf,'Renderer',rendertype);
axis equal
axis tight
box on
grid on

set(hFig1,'Units','Inches');
% set(hFig1, 'Position',fig_size);

xlabel('Global X-Dir', 'FontSize', 15);
ylabel('Global Y-Dir', 'FontSize', 15);
zlabel('Global Z-Dir', 'FontSize', 15);

end

