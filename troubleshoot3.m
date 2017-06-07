clc
clear



matPOINTS(:,:,1) = [0 2.5 0; 0 -2.5, 0];
matPOINTS(:,:,2) = [0 1.25 0; 0 -1.25 0];
matPOINTS(:,:,3) = [0.5 2.5 0; 0.5 -2.5 0];

% matPOINTS(:,:,1) = [0.5 2.5 0; 0 -2.5, 0];
% matPOINTS(:,:,2) = [0 1.25 0; 0 -1.25 0];
% matPOINTS(:,:,3) = [0 2.5 0; 0.5 -2.5 0];

% matPOINTS(:,:,1) = [0 1.25 0; 0 -2.5, 0];
% matPOINTS(:,:,2) = [0 2.5 0; 0 -1.25 0];
% matPOINTS(:,:,3) = [0.5 2.5 0; 0.5 -2.5 0];

% matPOINTS(:,:,1) = [0 0 0; 0 1 0];
% matPOINTS(:,:,2) = [1 1 0; 1 1 0];
% matPOINTS(:,:,3) = [1 0 0; 0 0 0];

matCOEFF = repmat([0 0 0 0 1 0],2,1);
vecUINF = [0.707 0 0.707];


[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matELOC, matPLEX, matDVECT, matALIGN, matVATT, matVNORM, matCENTER, matROTANG, matVSCOMB] = fcnTRIANG([], matPOINTS);

[hFig1] = fcnPLOTBODY(1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, vecUINF, matROTANG, [100 100 4 4], 'opengl');

% fpg = [0 0 1; 0.5 0 0.5];
fpg = repmat(matCENTER(2,:),2,1);
len = length(fpg(:,1));

dve_num = 2;
dir = 1;

for i = 1:3  
    edge_points = matVLST(matELST(matEIDX(dve_num,i),:),:);
    hold on
    if matVSCOMB(dve_num,i,dir) == 1
        plot3(edge_points(:,1), edge_points(:,2), edge_points(:,3),'--r','LineWidth',2)
    elseif matVSCOMB(dve_num,i,dir) == -1
        plot3(edge_points(:,1), edge_points(:,2), edge_points(:,3),'--b','LineWidth',2)
    end
    hold off
end

% set(hFig1,'Units','Inches');
% pos = get(hFig1,'Position');
% set(hFig1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(hFig1,'C:\Users\travis\Google Drive\PhD Proposal\figures\2dvesystem','-dpdf','-r0')

[q_ind] = fcnINDVEL([1:2], fpg, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, ones(len,1), matROTANG, matVSCOMB)

[q_ind1] = fcnINDVEL(ones(len,1), fpg, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, ones(len,1), matROTANG, matVSCOMB);
[q_ind2] = fcnINDVEL(ones(len,1)+1, fpg, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, ones(len,1), matROTANG, matVSCOMB);

hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind1(:,1), q_ind1(:,2), q_ind1(:,3), 1000, 'r')
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind2(:,1), q_ind2(:,2), q_ind2(:,3), 1000, 'b')
hold off

















