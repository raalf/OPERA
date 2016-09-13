% clc
clear

load('quad.mat');

fp = [-2 -2 0];

dvenum = [1]';

fpg = repmat(fp, length(dvenum),1,1);

COEFF = [0 1 0 1 0];

% dvenum = reshape(repmat(dvenum', length(fp), 1, 1),[],1,1);

COEFF = repmat(permute(COEFF(dvenum,:),[3 2 1]),3,1,1);

[q_ind] = fcnINDVEL(dvenum, fpg, COEFF, DVE, DVECT, VLST, DNORM, PLEX)
% q_ind = sum(q_ind,1)