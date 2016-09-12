% clc
clear

endpoints(:,:,1) = [0 0 0];
endpoints(:,:,2) = [1 1 0];
phi = pi/4;
yaw = pi/2;
fpl = [-6 1 0.5];
k = 1;

[al, bl, cl] = fcnVSIND(endpoints, phi, yaw, fpl, k);

%%

fpl = [-6 1 0.5];
phi = rad2deg(0.7853981);

orig = mean(endpoints,3);

k = 1;

Temp.DBL_EPS = 1e-14;
nuj = 0;
etaj = 0.5;

% assuming global ref frame is equal to local 
[a2x, b2x, c2x] = fcnVortexSheetInduction(Temp, fpl, orig, nuj, phi, etaj, k);