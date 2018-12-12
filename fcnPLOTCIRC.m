function [hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, colour, ppa)

for i = 1:valNELE
    corners = fcnGLOBSTAR(matVLST(matDVE(i,:),:) - matCENTER(i,:), [repmat(matROTANG(i,1),3,1) repmat(matROTANG(i,2),3,1) repmat(matROTANG(i,3),3,1)]);
    points = fcnPOLYGRID(corners(:,1), corners(:,2), ppa);
    
    len = size(points,1);

%     % points(:,2) is eta in local, points(:,1) is xsi
    circ = sum([0.5.*points(:,2).^2 points(:,2) 0.5.*points(:,1).^2 points(:,1) ones(size(points(:,1)))].*matCOEFF(i,:),2);
    vort = [matCOEFF(i,3).*points(:,1) + matCOEFF(i,4), matCOEFF(i,1).*points(:,2) + matCOEFF(i,2), points(:,2).*0];
    
    len = size(circ,1);
    tri = delaunay(points(:,1), points(:,2));
    
    circ_glob = fcnSTARGLOB([points circ], [repmat(matROTANG(i,1),len,1), repmat(matROTANG(i,2),len,1), repmat(matROTANG(i,3),len,1)]);
    circ_glob = circ_glob + matCENTER(i,:);
    
    vort_glob = fcnSTARGLOB(vort, [repmat(matROTANG(i,1),len,1) repmat(matROTANG(i,2),len,1) repmat(matROTANG(i,3),len,1)]);
    points_glob = fcnSTARGLOB([points points(:,1).*0], [repmat(matROTANG(i,1),len,1), repmat(matROTANG(i,2),len,1), repmat(matROTANG(i,3),len,1)]);
    points_glob = points_glob + matCENTER(i,:);
    
    hold on
%     trisurf(tri, circ_glob(:,1), circ_glob(:,2), circ_glob(:,3),'edgealpha',0,'facealpha',0.8);
    quiver3(points_glob(:,1), points_glob(:,2), points_glob(:,3), vort_glob(:,1), vort_glob(:,2), vort_glob(:,3))
    hold off
end





