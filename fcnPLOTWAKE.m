function [hFig1] = fcnPLOTWAKE(verbose, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, valWSIZE, matWVGRID, vecWDVESURFACE)

set(0,'CurrentFigure',hFig1);

col = single(vecWDVESURFACE)./max(single(vecWDVESURFACE));
hold on
colormap parula

% for j = 1:size(matWVGRID,2)
%     plot3(matWVLST(matWVGRID(1:end-valPRESTEPS,j),1),matWVLST(matWVGRID(1:end-valPRESTEPS,j),2),matWVLST(matWVGRID(1:end-valPRESTEPS,j),3),'-b','LineWidth',2);
% end

plotDVE = [1:valWNELE]';
patch('Faces',matWDVE(plotDVE,:,1),'Vertices',matWVLST,'FaceVertexCData',col(plotDVE),'FaceColor','flat','EdgeAlpha',0.4,'FaceAlpha',0.4);

if verbose == 1
    for ii = 1:valWNELE
        str = sprintf('%d',ii);
        text(matWCENTER(ii,1),matWCENTER(ii,2),matWCENTER(ii,3),str,'Color','k','FontSize',20);
    end
    
%     for ii = 1:length(matWVLST(:,1))
%         str = sprintf('%d',ii);
%         text(matWVLST(ii,1),matWVLST(ii,2),matWVLST(ii,3),str,'Color','g','FontSize',20);
%     end
%     
%     edge1 = matWVLST(matWELST(:,1),:);
%     edge2 = matWVLST(matWELST(:,2),:);
%     mid = (edge1+edge2)./2;
%     for ii = 1:length(mid)
%         str = sprintf('%d',ii);
%         text(mid(ii,1),mid(ii,2),mid(ii,3),str,'Color','b','FontSize',20);
%     end
%     
%     quiver3(matWCENTER(:,1),matWCENTER(:,2),matWCENTER(:,3), matWDVECT(:,1,1), matWDVECT(:,2,1), matWDVECT(:,3,1), 0.25, 'k'); % xsi
%     quiver3(matWCENTER(:,1),matWCENTER(:,2),matWCENTER(:,3), matWDVECT(:,1,2), matWDVECT(:,2,2), matWDVECT(:,3,2), 0.25, 'b'); % eta
%     quiver3(matWCENTER(:,1),matWCENTER(:,2),matWCENTER(:,3), matWDVECT(:,1,3), matWDVECT(:,2,3), matWDVECT(:,3,3), 0.25, 'm') % zeta (normal)
%     
end


hold off

end

