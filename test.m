clc
clear

load('matlab.mat');

hFig1 = figure(1);
clf(1);
patch('Faces',[1 2 3],'Vertices', matPLEX, 'FaceAlpha',0, 'linewidth',2)
grid minor
box on
axis tight
axis equal
xlabel('Local xsi','FontSize',15)
ylabel('Local eta','FontSize',15)
zlabel('Local zeta','FontSize',15)
for i = 1:3
    str = ['P' num2str(i)];
    text(matPLEX(i,1), matPLEX(i,2), matPLEX(i,3),str,'Color','b','FontSize',20);
end

p1 = matPLEX(1,:);
p2 = matPLEX(2,:);
p3 = matPLEX(3,:);
% Edge normal
N = nan(1,3,3);
N(:,:,1) = cross([0 0 1], p2-p1, 2);
N(:,:,2) = cross([0 0 1], p3-p2, 2);
N(:,:,3) = cross([0 0 1], p1-p3, 2);
N = N./sum(sqrt(N.^2),2);

comb = dot(N, repmat([0 1 0],1,1,3), 2);
comb = comb./abs(comb);
comb(isnan(comb)) = 0;
comb