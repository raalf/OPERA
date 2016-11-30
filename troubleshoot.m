% clc
clear

dvenum = 1;
fpg = [1 1 1];
% fpg =  [-0.207106781186548   0.207106781186548                   0];

matVLST = [...
            0.5 -0.5 0; ...
            0.5 0.5 0; ...
            -0.5 -0.5 0; ...
            -0.5 0.5 0 ...
            ];

%%
matDVE = [2 3 4];

matDVECT(:,:,1) = [-1 0 0];
matDVECT(:,:,2) = [0 -1 0];
matDVECT(:,:,3) = [0 0 1];

matPLEX = [0 0 0; ...
            1.000000000000000   1.000000000000000                   0; ...
            1.000000000000000                   0                   0];

[a1, a2, b1, b2, c3] = fcnHDVEIND2(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX)     

%%
matDVE = [2 4 3];

matDVECT(:,:,1) = [-0.707106781186547  -0.707106781186547                   0];
matDVECT(:,:,2) = [-0.707106781186547   0.707106781186547                   0];
matDVECT(:,:,3) = [0 0 -1];

matPLEX = [0 0 0; ...
            0.707106781186548   0.707106781186547                   0; ...
            1.414213562373095                   0                   0];

[a1, a2, b1, b2, c3] = fcnHDVEIND2(dvenum, fpg, matDVE, matDVECT, matVLST, matPLEX)