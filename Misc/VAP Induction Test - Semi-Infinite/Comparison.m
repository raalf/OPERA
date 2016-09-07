% clc
clear

endpoints(:,:,1) = [0 0 0];
endpoints(:,:,2) = [0 1 0];

fpg = [3 2 2];

orig = mean(endpoints,3);

k = 0;

%%
Temp.DBL_EPS = 1e-14;
nuj = 0;
phij = 0;
etaj = 0.5;

[a2x, b2x, c2x] = fcnVortexSheetInduction(Temp, fpg, orig, nuj, phij, etaj, k)

%%

% assuming global ref frame is equal to local 

[al, bl, cl] = fcnVSIND(endpoints, 0, [0 1 0], fpg, k)