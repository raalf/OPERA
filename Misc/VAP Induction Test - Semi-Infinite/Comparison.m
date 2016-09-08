clc
clear

endpoints(:,:,1) = [0 0 0];
endpoints(:,:,2) = [0 1 0];

fpg = [0 0 0; 1 2 5; 0.5 0.5 0; 0.5 1 0];
phi = [0 0 0 0]';


endpoints = repmat(endpoints,4,1,1);

orig = mean(endpoints,3);

k = [10 10 10 10]';
eta_VS = [0 1 0; 0 1 0; 0 1 0; 0 1 0];
% eta_VS = [1 0 0; 1 0 0; 1 0 0; 1 0 0];

[al, bl, cl] = fcnVSIND(endpoints, phi, eta_VS, fpg, k);

Temp.DBL_EPS = 1e-14;
nuj = 0;
phij = 0;
etaj = 0.5;

for i = 1:4
[a2x(i,:), b2x(i,:), c2x(i,:)] = fcnVortexSheetInduction(Temp, fpg(i,:), orig(i,:), nuj, phij, etaj, k(i));
end

al
a2x

bl
b2x

cl
c2x

%%

% clc
clear

endpoints(:,:,1) = [0 0 0];
endpoints(:,:,2) = [1 0 0];

fpg = [0 0 0; 1 2 5; 0.5 0.5 0; 0.5 1 0];
phi = [0 0 0 0]';


endpoints = repmat(endpoints,4,1,1);

orig = mean(endpoints,3);

k = [10 10 10 10]';
% eta_VS = [0 1 0; 0 1 0; 0 1 0; 0 1 0];
eta_VS = [1 0 0; 1 0 0; 1 0 0; 1 0 0];

Temp.DBL_EPS = 1e-14;
nuj = 0;
phij = 0;
etaj = 0.5;

% assuming global ref frame is equal to local 
num = 2;
[a2x, b2x, c2x] = fcnVortexSheetInduction(Temp, fpg(num,:), orig(num,:), nuj, phij, etaj, k(num));
[al, bl, cl] = fcnVSIND(endpoints(num,:,:), phi(num), eta_VS(num,:), fpg(num,:), k(num));



