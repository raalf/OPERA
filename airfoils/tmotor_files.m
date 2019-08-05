clc
clear

r = 0.4572/2;
y = [0.034290 0.045720 0.057150 0.068580 0.080010 0.091440 0.102870 0.114300 0.125730 0.137160 0.148590 0.160020 0.171450 0.182880 0.194310 0.205740 0.217170 0.228600]';
len = length(y);

pos = [0.25 0.5 0.75];
data(:,:,1) = dlmread('TMotor_Airfoil_25percent.dat','',1,0);
data(:,:,2) = dlmread('TMotor_Airfoil_50percent.dat','',1,0);
data(:,:,3) = dlmread('TMotor_Airfoil_75percent.dat','',1,0);

data(:,3,:) = repmat(reshape(pos,1,1,3), size(data(:,:,1),1),1);

yr = y./r;

hFig1 = figure(1);
clf(1);

plot3(data(:,1,1), data(:,2,1), data(:,3,1),'-k');
hold on
plot3(data(:,1,2), data(:,2,2), data(:,3,2),'-b');
plot3(data(:,1,3), data(:,2,3), data(:,3,3),'-m');
hold off

grid minor
box on
axis equal
axis tight

for j = 1:len
    x = permute(data(:,1:2,:),[3 1 2]);
    foil(:,:,j) = [permute(interp1(pos, x, yr(j), 'linear','extrap'), [2 3 1]) repmat(yr(j), size(data(:,:,1),1),1)];
    hold on
    plot3(foil(:,1,j), foil(:,2,j), foil(:,3,j),'-.r');
    hold off
    dlmwrite(['TM_S', num2str(j), '_', num2str(len), '.dat'], foil(:,1:2,j), '\t')
end



% out = reshape(permute(data, [2 1 3]), size(data, 2), [])';
% X = out(:,1);
% Y = out(:,2);
% Z = out(:,3);