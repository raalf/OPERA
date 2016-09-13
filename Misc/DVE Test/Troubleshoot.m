% % clc
clear

load('quad.mat');

fp = [.5 .5 .5];

dvenum = [1 2]';

fpg = repmat(fp, length(dvenum),1,1);

COEFF = [0 1 0 0 0; 0 1 0 0 0];

% dvenum = reshape(repmat(dvenum', length(fp), 1, 1),[],1,1);

COEFF = repmat(permute(COEFF(dvenum,:),[3 2 1]),3,1,1);

[q_ind] = fcnINDVEL(dvenum, fpg, COEFF, DVE, DVECT, VLST, DNORM, PLEX)
% q_ind = sum(q_ind,1)

%%
% clear
% 
% load('quad.mat');
% 
% fp = [0 -1 0];
% 
% dvenum = [1]';
% 
% fpg = repmat(fp, length(dvenum),1,1);
% 
% COEFF = [0 1 0 0 0];
% 
% % dvenum = reshape(repmat(dvenum', length(fp), 1, 1),[],1,1);
% 
% COEFF = repmat(permute(COEFF(dvenum,:),[3 2 1]),3,1,1);
% 
% [q_ind] = fcnINDVEL(dvenum, fpg, COEFF, DVE, DVECT, VLST, DNORM, PLEX)
% % q_ind = sum(q_ind,1)