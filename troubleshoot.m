clc
clear

dvenum = 1;
% fpg = [1 1 1];
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

x = -1:granularity:2;
y = -1:granularity:2;
z = -1:granularity:1;

[X,Y,Z] = meshgrid(x,y,z);

fpg = [reshape(X, [], 1, 1) reshape(Y, [], 1, 1) reshape(Z, [], 1, 1)];
fpg = unique(fpg,'rows');

% fpg = [4 2 1; 2 10 6];

len = length(fpg(:,1));

dvenum = ones(len,1);

%%
matDVE = [2 3 4];

matDVECT(:,:,1) = [1 0 0];
matDVECT(:,:,2) = [0 1 0];
matDVECT(:,:,3) = [0 0 1];

% ROLL = acos(dot(matDVECT(:,:,1),repmat([1 0 0],1,1,1),2));
% PITCH = acos(dot(matDVECT(:,:,2),repmat([0 1 0],1,1,1),2));
% YAW = acos(dot(matDVECT(:,:,3),repmat([0 0 1],1,1,1),2));

ROLL = -atan2(matDVECT(:,2,3), matDVECT(:,3,3));
PITCH = asin(matDVECT(:,1,3));
YAW = acos(dot(matDVECT(:,:,2),repmat([0 1 0],1,1,1),2));

matROTANG(:,1) = ROLL;
matROTANG(:,2) = PITCH;
matROTANG(:,3) = YAW;


matPLEX = [0 0 0; ...
            1.000000000000000   1.000000000000000                   0; ...
            1.000000000000000                   0                   0];

% [a1, a2, b1, b2, c3] = fcnHDVEIND2(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX);
[a1, a2, b1, b2, c3] = fcnHDVEIND(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX, matROTANG);

matCOEFF = [0 1 0 1 0];
D = [a1 a2 b1 b2 c3];
D = reshape(reshape(D', 1, 15, []), 3, 5, len);
q_ind = permute(sum(D.*repmat(reshape(matCOEFF(dvenum,:)',1,5,[]),3,1,1),2),[2 1 3]);
q_ind = reshape(permute(q_ind,[3 1 2]),[],3,1)./(-4*pi);



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
        

ROLL = -atan2(matDVECT(:,2,3), matDVECT(:,3,3));
PITCH = asin(matDVECT(:,1,3));
YAW = acos(dot(matDVECT(:,:,2),repmat([0 1 0],1,1,1),2));

matROTANG(:,1) = ROLL;
matROTANG(:,2) = PITCH;
matROTANG(:,3) = YAW;


% [a1, a2, b1, b2, c3] = fcnHDVEIND2(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX);
[a1, a2, b1, b2, c3] = fcnHDVEIND(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX, matROTANG);

matCOEFF = [0 0 0 1 0];
D = [a1 a2 b1 b2 c3];
D = reshape(reshape(D', 1, 15, []), 3, 5, len);
q_ind = permute(sum(D.*repmat(reshape(matCOEFF(dvenum,:)',1,5,[]),3,1,1),2),[2 1 3]);
q_ind = reshape(permute(q_ind,[3 1 2]),[],3,1)./(-4*pi);



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

%%
matDVE = [4 2 3];

matDVECT(:,:,1) = [0 1 0];
matDVECT(:,:,2) = [-1 0 0];
matDVECT(:,:,3) = [0 0 1];

ROLL = -atan2(matDVECT(:,2,3), matDVECT(:,3,3));
PITCH = asin(matDVECT(:,1,3));
YAW = acos(dot(matDVECT(:,:,2),repmat([0 1 0],1,1,1),2));

matROTANG(:,1) = ROLL;
matROTANG(:,2) = PITCH;
matROTANG(:,3) = YAW;

matPLEX = [0 0 0; ...
            0   1.000000000000000                   0; ...
            1.000000000000000                   0                   0];

% [a1, a2, b1, b2, c3] = fcnHDVEIND2(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX);
[a1, a2, b1, b2, c3] = fcnHDVEIND(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX, matROTANG);

matCOEFF = [0 -1 0 1 0];
D = [a1 a2 b1 b2 c3];
D = reshape(reshape(D', 1, 15, []), 3, 5, len);
q_ind = permute(sum(D.*repmat(reshape(matCOEFF(dvenum,:)',1,5,[]),3,1,1),2),[2 1 3]);
q_ind = reshape(permute(q_ind,[3 1 2]),[],3,1)./(-4*pi);


hFig7 = figure(7);
clf(7);

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