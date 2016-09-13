clc
clear

load('trouble2.mat')

% [~, bl, cl] = fcnVSIND2(endpoints, phi, yaw, fpl, k);

% bl(3,:)
% cl(3,:)

[~, bl2, cl2] = fcnVSIND2(endpoints(3,:,:), phi(3), yaw(3), fpl(3,:), k(3))

%%

endpoints3(:,:,1) = [0 0 0];
endpoints3(:,:,2) = [1 1 0];
phi3 = pi/4;
yaw3 = 0;
fpl3 = [.5 .5 .5];
k3 = 1;

% fpl3 = fpl(3,:);
% endpoints3 = endpoints(3,:,:)

% [al, bl, cl] = fcnVSIND(endpoints, phi, yaw, fpl, k)
[~, bl3, cl3] = fcnVSIND2(endpoints3, phi3, yaw3, fpl3, k3)