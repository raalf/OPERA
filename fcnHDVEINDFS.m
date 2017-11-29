function [infl_glob] = fcnHDVEINDFS(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX, dvetype, matROTANG, matVSCOMB, matCENTER)

flagPLOT = 1;

dvenum = reshape(dvenum, [], 1, 1); % Ensuring dvenum is a column vector

len = length(dvenum);

pa = fcnGLOBSTAR(fpg - matCENTER(dvenum,:), matROTANG(dvenum,1), matROTANG(dvenum,2), matROTANG(dvenum,3));

p1 = permute(matPLEX(1,:,dvenum),[3 2 1]);
p2 = permute(matPLEX(2,:,dvenum),[3 2 1]);
p3 = permute(matPLEX(3,:,dvenum),[3 2 1]);

fn = cross(p3 - p1, p2 - p1);

N = zeros(len,3,3);
N(:,:,1) = cross(fn, p2-p1, 2);
N(:,:,2) = cross(fn, p3-p2, 2);
N(:,:,3) = cross(fn, p1-p3, 2);
N = N./sum(sqrt(N.^2),2);

%% (b1,b2,b3)
% b1 = (p3 - p1);
% b3 = cross(p3 - p1, p2 - p1);
% b2 = cross(b3, b1);
% b1 = b1./sqrt(sum(b1.^2,2));
% b2 = b2./sqrt(sum(b2.^2,2));
% b3 = b3./sqrt(sum(b3.^2,2));

b1 = matDVECT(dvenum,:,1);
b2 = matDVECT(dvenum,:,2);
b3 = matDVECT(dvenum,:,3);

%% Projection of P_A on plane
pb = pa - dot((pa - p1), b3, 2).*b3;

z = dot((pa - pb), b3, 2);

%% Influence from S1 S2 and S3
q1 = p1 - pb;
q2 = p2 - pb;
q3 = p3 - pb;

c3 = b3;

delta = 0.2;
h = sqrt(z.^2 + delta.^2);

b(:,:,1) = b1;
b(:,:,2) = b2;
b(:,:,3) = b3;

if flagPLOT == 1 && len == 1
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
    
    h1 = quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), matDVECT(:,1,1), matDVECT(:,2,1), matDVECT(:,3,1), 0.25, 'k'); % xsi
    h2 = quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), matDVECT(:,1,2), matDVECT(:,2,2), matDVECT(:,3,2), 0.25, 'b'); % eta
    h3 = quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), matDVECT(:,1,3), matDVECT(:,2,3), matDVECT(:,3,3), 0.25, 'm'); % zeta (normal)
    
    scatter3(fpg(:,1), fpg(:,2), fpg(:,3), 100, 'ok', 'filled')
    text(fpg(:,1), fpg(:,2), fpg(:,3), 'PA', 'FontSize',20,'Color','m')
    
    hold off
    axis equal
    box on
    grid minor
    
    title('Global Coordinate System','FontSize',15);
    xlabel('Global X-Dir', 'FontSize', 15);
    ylabel('Global Y-Dir', 'FontSize', 15);
    zlabel('Global Z-Dir', 'FontSize', 15);
    
    % LOCAL ------------------------------------------------------------------------------------------------------------------------------
    ax(2) = subplot(1,2,2);
    patch(matPLEX(:,1), matPLEX(:,2), matPLEX(:,3),'LineWidth',2,'EdgeColor','k','FaceAlpha',0);
    
    hold on
    scatter3(pa(:,1), pa(:,2), pa(:,3), 100, 'ok', 'filled')
    text(pa(:,1), pa(:,2), pa(:,3), 'PA', 'FontSize',20,'Color','m')
    
    scatter3(pb(:,1), pb(:,2), pb(:,3), 100, 'og', 'filled')
    text(pb(:,1), pb(:,2), pb(:,3), 'PB', 'FontSize',20,'Color','m')

    % (b1,b2,b3)
    quiver3(pb(:,1),pb(:,2),pb(:,3), b1(1,1), b1(1,2), b1(1,3), 0.2, 'b', 'LineWidth',2)
    quiver3(pb(:,1),pb(:,2),pb(:,3), b2(1,1), b2(1,2), b2(1,3), 0.2, 'b', 'LineWidth',2)
    quiver3(pb(:,1),pb(:,2),pb(:,3), b3(1,1), b3(1,2), b3(1,3), 0.2, 'b', 'LineWidth',2)

    for i = 1:3
        text(matPLEX(i,1), matPLEX(i,2), matPLEX(i,3), ['P' num2str(i)], 'FontSize',20,'Color','r')
    end
    
    hold off
    axis equal
    box on
    grid minor
    
    title('Local Coordinate System','FontSize',15);
    xlabel('Local Xsi-Dir', 'FontSize', 15);
    ylabel('Local Eta-Dir', 'FontSize', 15);
    zlabel('Local Zeta-Dir', 'FontSize', 15);
    
    view([0 90])
        
%     linkprop(ax,{'xlim','ylim','zlim','CameraUpVector', 'CameraPosition', 'CameraTarget'});   
    
end

% S1 (pb,p1,p2)
[infl_1] = fcnSNINF(q1, q2, c3, b, z, h);

% S2 (pb,p2,p3)
[infl_2] = fcnSNINF(q2, q3, c3, b, z, h);

% S3 (pb,p3,p1)
[infl_3] = fcnSNINF(q3, q1, c3, b, z, h);

%% Combining S1, S2, S3 to make S

comb = [dot(N(:,:,1), q1, 2) dot(N(:,:,2),q2,2) dot(N(:,:,3),q3,2)];
comb = comb./abs(comb);
comb = reshape(comb',1,3,[]);
comb(isnan(comb)) = 0;

infl_1 = infl_1.*repmat(permute(comb(:,1,:),[2 1 3]),3,6,1);
infl_2 = infl_2.*repmat(permute(comb(:,2,:),[2 1 3]),3,6,1);
infl_3 = infl_3.*repmat(permute(comb(:,3,:),[2 1 3]),3,6,1);

infl = infl_1 + infl_2 + infl_3;

dvenum = reshape(repmat(dvenum,1,6,1)',[],1,1);

infl_glob = fcnSTARGLOB(reshape(permute(infl,[2 3 1]),[],3,1), matROTANG(dvenum,1), matROTANG(dvenum,2), matROTANG(dvenum,3));
infl_glob = reshape(infl_glob',3,6,[]);

infl_glob(isnan(infl_glob)) = 0;

end
