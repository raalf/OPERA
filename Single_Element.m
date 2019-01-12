clear
% clc

%% Preamble

strFILE = 'inputs/2dve.dat';

% [~, strATYPE, vecSYM, ~, ~, ~, valALPHA, valBETA, ~, ~] = fcnOPREAD(strFILE);
strATYPE = 'what';
valALPHA = 0;
valBETA = 0;
vecSYM = [];


% matPOINTS(:,:,1) = [0 0 0];
% matPOINTS(:,:,2) = [0  0.5 0];
% matPOINTS(:,:,3) = [0.5  0.5 0];

matPOINTS(:,:,1) = [0.5 0 0];
matPOINTS(:,:,2) = [0  0 0];
matPOINTS(:,:,3) = [0  0.5 0];


[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
    matELOC, matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG, matCONTROL] = fcnTRIANG(matPOINTS);

vecUINF = fcnUINFWING(valALPHA, 0);

%% fpl
fpg = [4 5 0]
fpl = fcnGLOBSTAR(fpg - matCENTER, matROTANG);

%% Coefficients
matCOEFF = [1 1 1 1 1];

syms xi eta real
mu = 0.5.*matCOEFF(1).*eta.^2 + matCOEFF(2).*eta + 0.5.*matCOEFF(3).*xi.^2 + matCOEFF(4).*xi + matCOEFF(5);

s = [xi eta xi.*0];
r = fpl;
n = matDVECT(:,:,3);

term = repmat(mu,1,3).*((n./(sqrt(sum((r - s).^2,2)).^3)) - dot(3.*repmat(n,length(xi(:)),1), r - s, 2).*((r - s)./(sqrt(sum((r - s).^2,2)).^5)));

xi1 = matPLEX(1,1);
xi2 = matPLEX(2,1);
xi3 = matPLEX(3,1);
eta1 = matPLEX(1,2);
eta2 = matPLEX(2,2);
eta3 = matPLEX(3,2);

le = (((eta3 - eta2).*xi)./(xi3 - xi2)) + eta2 - ((xi2.*(eta3 - eta2))./(xi3 - xi2));
te = (((eta3 - eta1).*xi)./(xi3 - xi1)) + eta1 - ((xi1.*(eta3 - eta1))./(xi3 - xi1));

nint = vpaintegral(int(term, eta, le, te), xi, xi1, xi3);
% v = fcnSTARGLOB((1./(4*pi)).*nint, matROTANG)
infl_loc = fcnHDVEIND_DB(1, 1, fpg, matPLEX, matROTANG, matCENTER);
q_ind = permute(sum(infl_loc.*repmat(reshape(matCOEFF(1,:)',1,5,[]),3,1,1),2),[2 1 3]);

abs((nint - q_ind)./q_ind).*100
% [q_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCONTROL)

%% Plot

% vecUINF = [cosd(10) 0 sind(10)]

% [hFig1] = fcnPLOTBODY(1, matDVE, valNELE, matVLST, matELST, matDVECT, matCONTROL, matPLEX, [], vecUINF, matROTANG, [3 1 4 4], 'opengl');
% % [hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), vecUINF, matROTANG, 'r', 10);
% % view([-30 17])
% 
% % granularity = 0.05;
% % y = [-0.2:granularity:0.8];
% % x = [0.05, matCENTER(:,1), 0.4];
% % z = [-0.2:granularity:0.2];
% 
% % granularity = 0.025;
% % x = [-0.2:granularity:0.8];
% % y = [0.2, matCENTER(:,2), 0.45];
% % z = [-0.2:granularity:0.2];
% 
% 
% 
% % granularity = 0.025;
% % x = [-0.1:granularity:0.6];
% % y = [-0.2:granularity:0.6];
% % z = x.*0;
% 
% % granularity = 0.01;
% % z = -.5:granularity:.5;
% % len = length(z);
% % y = zeros(len,1) + matCENTER(:,2);
% % x = zeros(len,1) + matCENTER(:,1);
% 
% granularity = 0.01;
% z = -.5:granularity:.5;
% len = length(z);
% % y = zeros(len,1) + matCENTER(:,2);
% y = 0.01:0.1:0.49;
% x = -0.05:0.06:0.35;
% 
% 
% % z(z==0) = [];
% [X,Y,Z] = meshgrid(x,y,z);
% fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');
% 
% % fpg = [matCENTER; matCENTER + [0.1 0 0]];
% % fpg = [matCENTER; matCENTER + [0 0 0.1]];
% % fpg = matCENTER + [2 1 1];
% 
% % fpg = [0.6 0.5 1.5; 0.23 0.5 0.6]
% % fpg = [-0.1000   -0.1000   -0.2550]
% % fpg = [0    0.1667   -0.2000]
% % fpg = matCENTER
% 
% [q_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCONTROL);
% 
% figure(1);
% hold on
% quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3),'b')
% hold off







