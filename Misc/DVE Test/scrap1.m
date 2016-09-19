clc
clear

v1 = [0.5 2 5];

DVECT(:,:,1) = [0 1 0];
DVECT(:,:,2) = [1 0 0];
DVECT(:,:,3) = [0 0 -1];

A = [DVECT(:,:,1)' DVECT(:,:,2)' DVECT(:,:,3)']

v2 = A*v1'