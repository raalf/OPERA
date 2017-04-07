clc
clear

dvenum = 1;
% fpg = [1 2 1];
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
% fpg = [-0.5 1 1]
% fpg = [0.5 0.5 1]

len = length(fpg(:,1));

dvenum = ones(len,1);

test = 1;

matCOEFF = [0 0 0 1 0];

%
matDVE = [2 3 4];

matDVECT(:,:,1) = [1 0 0];
matDVECT(:,:,2) = [0 1 0];
matDVECT(:,:,3) = [0 0 1];

ROLL = -atan2(matDVECT(:,2,3), matDVECT(:,3,3));
PITCH = asin(matDVECT(:,1,3));
YAW = acos(dot(matDVECT(:,:,2),repmat([0 1 0],1,1,1),2));

matROTANG(:,1) = ROLL;
matROTANG(:,2) = PITCH;
matROTANG(:,3) = YAW;


matPLEX = [0 0 0; ...
            1.000000000000000   1.000000000000000                   0; ...
            1.000000000000000                   0                   0];

        
[q_ind] = fcnSDVEVEL(fpg, 1, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG);

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