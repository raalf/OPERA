function [velocities_out, speeds_out, points_out] = fcnSURVEL(matDVE, valNELE, matVLST, matCENTER, matCOEFF, matUINF, matROTANG, ppa)

    velocities_out = [];
    speeds_out = [];
    points_out = [];

for i = 1:valNELE
    corners = fcnGLOBSTAR(matVLST(matDVE(i,:),:) - matCENTER(i,:), [repmat(matROTANG(i,1),3,1) repmat(matROTANG(i,2),3,1) repmat(matROTANG(i,3),3,1)]);
    points = fcnPOLYGRID(corners(:,1), corners(:,2), ppa);
    
%     % points(:,2) is eta in local, points(:,1) is xsi
    vort = [matCOEFF(i,3).*points(:,1) + matCOEFF(i,4) + matCOEFF(i,5).*points(:,2), matCOEFF(i,1).*points(:,2) + matCOEFF(i,2) + matCOEFF(i,5).*points(:,1), points(:,2).*0];
    len = size(vort,1);
    vort_glob = fcnSTARGLOB(vort, [repmat(matROTANG(i,1),len,1) repmat(matROTANG(i,2),len,1) repmat(matROTANG(i,3),len,1)]);

    velocities = (vort_glob.*2) + repmat(matUINF(i,:),len,1);
    speeds = sqrt(sum(velocities.^2,2));
    
    tri = delaunay(points(:,1), points(:,2));
    
    hold on
    trisurf(tri, points(:,1), points(:,2), points(:,1).*0, speeds, 'edgealpha',0,'facealpha',0.8);
%     quiver3(points_glob(:,1), points_glob(:,2), points_glob(:,3), vort_glob(:,1), vort_glob(:,2), vort_glob(:,3))
    hold off
    
%     velocities_out(:,:,i) = velocities;
%     speeds_out(:,1,i) = speeds;
%     points_out(:,:,i) = points;
    
    velocities = [];
    speeds = [];
    points = [];
    
end





