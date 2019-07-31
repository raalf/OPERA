clear
% clc

%% Preamble
strFILE = 'inputs/2dve.dat';
strATYPE = 'what';
valALPHA = 0;
valBETA = 0;
vecSYM = [];

matPOINTS(:,:,1) = [1 0 0];
matPOINTS(:,:,2) = [0  0 0];
% matPOINTS(:,:,3) = [0.9  0.5 0];
xp = 0
% xp = 0.99
% xp = 0.5

matPOINTS(:,:,3) = [xp 0.5 0];
[TR, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
    matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG, matCONTROL] = fcnTRIANG(matPOINTS, false);

% matPOINTS(:,:,3) = [xp  -0.5 0];
% [TR, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
%     matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG, matCONTROL] = fcnTRIANG(matPOINTS,true);

%%
vecDVESYM = false(size(matCENTER, 1), 1);
vecUINF = fcnUINFWING(valALPHA, 0);

% matCOEFF = [0 0 0   -0.3236    0.0000   -0.1021];
% matCOEFF = -[10 0 -20 1 0 1];
% matCOEFF = -[0 0 0 0 1 0];
matCOEFF = -[1 1 1 1 1 1];

%% Plot
granularity = 0.001;
x = matCENTER(1);
y = matCENTER(2);
% x = 0;
% y = 0.5;
% x = 0.25;
% y = 0.25;
x = 1;
y = 2;
z = [-2:granularity:2];

[X,Y,Z] = meshgrid(x,y,z);
fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');

% fpg = [x y 0.187]

%% YEH
len = length(fpg(:,1));
dvenum = reshape(repmat(1:valNELE,len,1),[],1);
dvetype = ones(size(dvenum));
fpg = repmat(fpg,valNELE,1);
tic
infl_loc = fcnHDVEIND_DB_YEH(dvenum, dvetype, fpg, matPLEX, matROTANG, matCENTER, [], 0);
toc
q_ind = permute(sum(infl_loc.*repmat(reshape(matCOEFF(dvenum,:)',1,6,[]),3,1,1),2),[2 1 3]);
q_ind = reshape(permute(q_ind,[3 1 2]),[],3,1)./(4*pi);
q_ind = fcnSTARGLOB(q_ind, matROTANG(dvenum,:));
q_ind1 = reshape(sum(permute(reshape(q_ind',3,[],valNELE),[3 1 2]),1),3,[],1)';

%% FTJ
len = length(fpg(:,1));
dvenum = reshape(repmat(1:valNELE,len,1),[],1);
dvetype = ones(size(dvenum));
fpg = repmat(fpg,valNELE,1);
tic
infl_loc = fcnHDVEIND_DB(dvenum, dvetype, fpg, matPLEX, matROTANG, matCENTER, [], 5e-3);
toc
q_ind = permute(sum(infl_loc.*repmat(reshape(matCOEFF(dvenum,:)',1,6,[]),3,1,1),2),[2 1 3]);
q_ind = reshape(permute(q_ind,[3 1 2]),[],3,1)./(4*pi);
q_ind = fcnSTARGLOB(q_ind, matROTANG(dvenum,:));
q_ind2 = reshape(sum(permute(reshape(q_ind',3,[],valNELE),[3 1 2]),1),3,[],1)';

%%
figure(2);
clf(2);
subplot(3,1,1)
plot(q_ind1(:,1), z, '-k')
hold on
plot(q_ind2(:,1), z, '--b')
hold off
grid minor
box on
axis tight

subplot(3,1,2)
plot(q_ind1(:,2), z, '-k')
hold on
plot(q_ind2(:,2), z, '--b')
hold off
grid minor
box on
axis tight

subplot(3,1,3)
plot(q_ind1(:,3), z, '-k')
hold on
plot(q_ind2(:,3), z, '--b')
hold off
grid minor
box on
axis tight





