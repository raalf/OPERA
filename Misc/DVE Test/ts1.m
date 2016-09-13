% clc
clear

endpoints(:,:,1) = [0 0 0];
endpoints(:,:,2) = [1 1 0];
phi = pi/4;
yaw = 0;
fpl = [.5 .5 .5];
k = 1;

% [al, bl, cl] = fcnVSIND(endpoints, phi, yaw, fpl, k)
[~, bl, cl] = fcnVSIND2(endpoints, phi, yaw, fpl, k)

%%

endpoints(:,:,1) = [0 0 0];
endpoints(:,:,2) = [1 1 0];

fpl = [2.5 2.5000 3];
phi = rad2deg(pi/4);

orig = mean(endpoints,3);

k = 1;

Temp.DBL_EPS = 1e-14;
nuj = 0;
etaj = 0.5;

% assuming global ref frame is equal to local 
[a2x, b2x, c2x] = fcnVortexSheetInduction(Temp, fpl, orig, nuj, phi, etaj, k);
disp(b2x)
disp(c2x)

%%
% clc
clear

endpoints(:,:,1) = [0 0 0];
endpoints(:,:,2) = [1 0 0];
phi = 0;
yaw = 0;
fpl = [4 0 3];
k = 1;

[al, bl, cl] = fcnVSIND2(endpoints, phi, yaw, fpl, k)

%%
% clc
clear

endpoints(:,:,1) = [0 0 0];
endpoints(:,:,2) = [0 1 0];
phi = 0;
yaw = pi/2;
fpl = [0 4 3];
k = 1;

[al, bl, cl] = fcnVSIND2(endpoints, phi, yaw, fpl, k)



