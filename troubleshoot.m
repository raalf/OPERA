clc
clear

dvenum = 1;
fpg = [0 1 0];
% fpg =  [-0.207106781186548   0.207106781186548                   0];

% matVLST = [...
%             0.5 -0.5 0; ...
%             0.5 0.5 0; ...
%             -0.5 -0.5 0; ...
%             -0.5 0.5 0 ...
%             ];

matVLST = [...
            0 -0.5 0; ...
            0 0 0; ...
            1 1 0; ...
            1 0 0 ...
            ];
        
granularity = 0.25;

x = -2:granularity:2;
y = -2:granularity:2;
z = -2:granularity:2;

[X,Y,Z] = meshgrid(x,y,z);

fpg = [reshape(X, [], 1, 1) reshape(Y, [], 1, 1) reshape(Z, [], 1, 1)];
fpg = unique(fpg,'rows');

len = length(fpg(:,1));

dvenum = ones(len,1);

%%
matDVE = [2 3 4];

matDVECT(:,:,1) = [1 0 0];
matDVECT(:,:,2) = [0 1 0];
matDVECT(:,:,3) = [0 0 1];

matPLEX = [0 0 0; ...
            1.000000000000000   1.000000000000000                   0; ...
            1.000000000000000                   0                   0];

[a1, a2, b1, b2, c3] = fcnHDVEIND2(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX);
% [a1, a2, b1, b2, c3] = fcnHDVEIND(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX);

q_ind = a1 + a2 + b1 + b2 + c3;
q_ind = a2 + b2 + c3;

hFig5 = figure(5);
clf(5);

patch(matVLST(matDVE(:),1), matVLST(matDVE(:),2), matVLST(matDVE(:),3),'r','FaceAlpha',0.5);
hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3));
hold off

grid on
axis tight
box on
axis equal

xlabel('X-Dir','FontSize',15);
ylabel('Y-Dir','FontSize',15);
zlabel('Z-Dir','FontSize',15);

%%
matDVE = [2 4 3];

matDVECT(:,:,1) = [0.707106781186547  0.707106781186547                   0];
matDVECT(:,:,2) = [0.707106781186547  -0.707106781186547                   0];
matDVECT(:,:,3) = [0 0 -1];

matPLEX = [0 0 0; ...
            0.707106781186548   0.707106781186547                   0; ...
            1.414213562373095                   0                   0];

[a1, a2, b1, b2, c3] = fcnHDVEIND2(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX);
% [a1, a2, b1, b2, c3] = fcnHDVEIND(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX);

q_ind = a1 + a2 + b1 + b2 + c3;
q_ind = a2 + b2 + c3;

hFig6 = figure(6);
clf(6);

patch(matVLST(matDVE(:),1), matVLST(matDVE(:),2), matVLST(matDVE(:),3),'r','FaceAlpha',0.5);
hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3));
hold off

grid on
axis tight
box on
axis equal

xlabel('X-Dir','FontSize',15);
ylabel('Y-Dir','FontSize',15);
zlabel('Z-Dir','FontSize',15);