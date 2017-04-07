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
% fpg = [0.5 0.25 0.2]
% fpg = [1 -0.5 1];

len = length(fpg(:,1));

dvenum = ones(len,1);

test = 1;

matCOEFF = [0 1 0 0 0];

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

        
[a1, a2, b1, b2, c3] = fcnHDVEIND(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX, 1, matROTANG);

D = [a1 a2 b1 b2 c3];
D = reshape(reshape(D', 1, 15, []), 3, 5, len);

coeff = matCOEFF;

q_ind = permute(sum(D.*repmat(reshape(coeff(dvenum,:)',1,5,[]),3,1,1),2),[2 1 3]);
q_ind = reshape(permute(q_ind,[3 1 2]),[],3,1)./(-4*pi);

% q_ind = fcnROTVECT(1, q_ind, matDVECT);
% q_ind = fcnTOGLOB(ones(length(q_ind),1),q_ind, matDVE, matDVECT, matVLST);

% q_ind = fcnSTARGLOB(q_ind, matROTANG(:,1), matROTANG(:,2), matROTANG(:,3));

q_ind(test,:)

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
matDVE = [3 4 2];

matDVECT(:,:,1) = [-0.707106781186547  -0.707106781186547                   0];
matDVECT(:,:,2) = [0.707106781186547  -0.707106781186547                   0];
matDVECT(:,:,3) = [0 0 1];

matPLEX = [0 0 0; ...
            0.707106781186548   0.707106781186547                   0; ...
            1.414213562373095                   0                   0];
        
% ROLL = acos(dot(matDVECT(:,:,1),repmat([1 0 0],1,1,1),2));
% PITCH = acos(dot(matDVECT(:,:,2),repmat([0 1 0],1,1,1),2));
% YAW = acos(dot(matDVECT(:,:,3),repmat([0 0 1],1,1,1),2)); 

ROLL = -atan2(matDVECT(:,2,3), matDVECT(:,3,3));
PITCH = asin(matDVECT(:,1,3));
% YAW = -acos(dot(matDVECT(:,:,2),repmat([0 1 0],1,1,1),2));
YAW = atan2(matDVECT(:,2,1), matDVECT(:,1,1));

matROTANG(:,1) = ROLL;
matROTANG(:,2) = PITCH;
matROTANG(:,3) = YAW;

[a1, a2, b1, b2, c3] = fcnHDVEIND(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX, 1, matROTANG);

D = [a1 a2 b1 b2 c3];
D = reshape(reshape(D', 1, 15, []), 3, 5, len);

coeff = matCOEFF;

q_ind = permute(sum(D.*repmat(reshape(coeff(dvenum,:)',1,5,[]),3,1,1),2),[2 1 3]);
q_ind = reshape(permute(q_ind,[3 1 2]),[],3,1)./(-4*pi);

% q_ind = fcnSTARGLOB(q_ind, matROTANG(:,1), matROTANG(:,2), matROTANG(:,3));

q_ind(test,:)

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
% YAW = acos(dot(matDVECT(:,:,2),repmat([0 1 0],1,1,1),2));
% YAW = atan2(-matDVECT(:,2,1), -matDVECT(:,1,1))
YAW = atan2(matDVECT(:,2,1), matDVECT(:,1,1));

matROTANG(:,1) = ROLL;
matROTANG(:,2) = PITCH;
matROTANG(:,3) = YAW;

matPLEX = [0 0 0; ...
            0   1.000000000000000                   0; ...
            1.000000000000000                   0                   0];

[a1, a2, b1, b2, c3] = fcnHDVEIND(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX, 1, matROTANG);

D = [a1 a2 b1 b2 c3];
D = reshape(reshape(D', 1, 15, []), 3, 5, len);

coeff = matCOEFF;

q_ind = permute(sum(D.*repmat(reshape(coeff(dvenum,:)',1,5,[]),3,1,1),2),[2 1 3]);
q_ind = reshape(permute(q_ind,[3 1 2]),[],3,1)./(-4*pi);

% q_ind = fcnSTARGLOB(q_ind, matROTANG(:,1), matROTANG(:,2), matROTANG(:,3));

q_ind(test,:)

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