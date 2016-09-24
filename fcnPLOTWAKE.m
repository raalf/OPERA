function [hFig1] = fcnPLOTWAKE(verbose, WDVE, WNELE, WVLST, WELST, WDVECT, WCENTER)
%FCNPLOTWAKE Summary of this function goes here
%   Detailed explanation goes here

hFig1 = figure(1);

hold on

patch('Faces',WDVE(:,:,1),'Vertices',WVLST,'FaceColor','b','LineWidth',2);
alpha(0.5);

if verbose == 1
    for ii = 1:WNELE
        str = sprintf('%d',ii);
        text(WCENTER(ii,1),WCENTER(ii,2),WCENTER(ii,3),str,'Color','k','FontSize',20);
    end
    
    for ii = 1:length(WVLST(:,1))
        str = sprintf('%d',ii);
        text(WVLST(ii,1),WVLST(ii,2),WVLST(ii,3),str,'Color','g','FontSize',20);
    end
    
    edge1 = WVLST(WELST(:,1),:);
    edge2 = WVLST(WELST(:,2),:);
    mid = (edge1+edge2)./2;
    for ii = 1:length(mid)
        str = sprintf('%d',ii);
        text(mid(ii,1),mid(ii,2),mid(ii,3),str,'Color','b','FontSize',20);
    end
    
    quiver3(WCENTER(:,1),WCENTER(:,2),WCENTER(:,3), WDVECT(:,1,1), WDVECT(:,2,1), WDVECT(:,3,1), 0.25, 'b') % eta
    quiver3(WCENTER(:,1),WCENTER(:,2),WCENTER(:,3), WDVECT(:,1,2), WDVECT(:,2,2), WDVECT(:,3,2), 0.25, 'k') % xi
    quiver3(WCENTER(:,1),WCENTER(:,2),WCENTER(:,3), WDVECT(:,1,3), WDVECT(:,2,3), WDVECT(:,3,3), 0.25, 'm') % zeta (normal)
    
end


hold off

end

