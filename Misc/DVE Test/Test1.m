clc
clear

% insert hdve num and fpg and get back induced velocity vector

load('quad.mat');

gran = 0.5;
[X Y Z] = meshgrid(-3:gran:3);

fp = [reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)];


% scatter3(fp(:,1), fp(:,2), fp(:,3),'ok')

%%

% fp(:,1) = 0;

% fp = unique(fp,'rows');

points = length(fp);

dvenum = [1 2]';

fpg = repmat(fp, length(dvenum),1,1);


COEFF = [0 1 0 1 0; 0 1 0 1 0];

dvenum = reshape(repmat(dvenum', length(fp), 1, 1),[],1,1);

COEFF = repmat(permute(COEFF(dvenum,:),[3 2 1]),3,1,1);



[q_ind] = fcnINDVEL(dvenum, fpg, COEFF, DVE, DVECT, VLST, DNORM, PLEX);

q_ind = q_ind(1:points,:) + q_ind(points+1:end,:);


[hFig1] = fcnPLOTBODY(1, DVE, NELE, VLST, ELST, DVECT);

hold on
quiver3(fp(:,1), fp(:,2), fp(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3),0)
hold off
