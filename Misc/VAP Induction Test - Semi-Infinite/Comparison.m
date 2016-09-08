clc
clear

endpoints(:,:,1) = [0 0 0];
endpoints(:,:,2) = [0 1 0];

fpg = [0 0 0; 1 2 5; 0.5 0.5 0; 0.5 1 0];
phi = [40 40 40 40]';


endpoints = repmat(endpoints,4,1,1);

orig = mean(endpoints,3);

k = [10 10 10 10]';
eta_VS = [0 1 0; 0 1 0; 0 1 0; 0 1 0];
% eta_VS = [1 0 0; 1 0 0; 1 0 0; 1 0 0];

[al, bl, cl] = fcnVSIND(endpoints, deg2rad(phi), eta_VS, fpg, k);

Temp.DBL_EPS = 1e-14;
nuj = 0;
etaj = 0.5;

for i = 1:4
[a2x(i,:), b2x(i,:), c2x(i,:)] = fcnVortexSheetInduction(Temp, fpg(i,:), orig(i,:), nuj, phi(i), etaj, k(i));
end

al
a2x

bl
b2x

cl
c2x

%%

clc
clear

endpoints(:,:,1) = [0 0 0];
endpoints(:,:,2) = [0 1 0];

fpg = [1 2 5];
phi = 0;

orig = mean(endpoints,3);

k = 10;
% eta_VS = [0 1 0; 0 1 0; 0 1 0; 0 1 0];
eta_VS = [0 1 0];

Temp.DBL_EPS = 1e-14;
nuj = 0;
etaj = 0.5;

% assuming global ref frame is equal to local 
[a2x, b2x, c2x] = fcnVortexSheetInduction(Temp, fpg, orig, nuj, phi, etaj, k);
[al, bl, cl] = fcnVSIND(endpoints, deg2rad(phi), eta_VS, fpg, k);

al
a2x

bl
b2x

cl
c2x


