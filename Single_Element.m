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
% fpg = [-0 0.2 0]

% Wide
% granularity = 0.05;
% x = [-5:granularity:7];
% y = [-5:granularity:5];
% z = [-5:granularity:5];
% [X,Y,Z] = meshgrid(x,y,z);
% fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');

% Narrow
granularity = 0.005;
x = [-0.05:granularity:0.55];
y = [-0.05:granularity:0.55];
z = [-0.1:granularity:0.1];
[X,Y,Z] = meshgrid(x,y,z);
fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');

% fpg = [0.1183333333333333597936487535662308800966 0.1416666666666666629659232512494781985879 0.0000000000000000000000000000000000000000]

fpl = fcnGLOBSTAR(fpg - matCENTER, matROTANG);
len = size(fpl,1);

%% Coefficients
matCOEFF = [1 1 1 1 1];

syms xi eta real
mu = 0.5.*matCOEFF(1).*eta.^2 + matCOEFF(2).*eta + 0.5.*matCOEFF(3).*xi.^2 + matCOEFF(4).*xi + matCOEFF(5);

s = [xi eta xi.*0];
n = matDVECT(:,:,3);

xi1 = matPLEX(1,1);
xi2 = matPLEX(2,1);
xi3 = matPLEX(3,1);
eta1 = matPLEX(1,2);
eta2 = matPLEX(2,2);
eta3 = matPLEX(3,2);

le = (((eta3 - eta2).*xi)./(xi3 - xi2)) + eta2 - ((xi2.*(eta3 - eta2))./(xi3 - xi2));
te = (((eta3 - eta1).*xi)./(xi3 - xi1)) + eta1 - ((xi1.*(eta3 - eta1))./(xi3 - xi1));

infl_loc = fcnHDVEIND_DB(ones(len,1), ones(len,1), fpg, matPLEX, matROTANG, matCENTER);
q_ind = permute(sum(infl_loc.*repmat(reshape(matCOEFF(1,:)',1,5,[]),3,1,1),2),[2 1 3]);

if len > 1
    delete('BadFPs.txt');
end

for i = 1:size(fpl,1)
% parfor i = 1:size(fpl,1)
    try
        r = fpl(i,:);
        term = repmat(mu,1,3).*((n./(sqrt(sum((r - s).^2,2)).^3)) - dot(3.*repmat(n,length(xi(:)),1), r - s, 2).*((r - s)./(sqrt(sum((r - s).^2,2)).^5)));
        nint = vpaintegral(int(term, eta, le, te), xi, xi1, xi3);
        
        d = abs((nint - q_ind(:,:,i))./q_ind(:,:,i)).*100;
        
        if len > 1
            if max(real(d)) > 1e-3
                disp([num2str(i),': Bad. FPL = [',num2str(fpl(i,1)),', ', num2str(fpl(i,2)), ', ', num2str(fpl(i,3)), ']'])
                dlmwrite('BadFPs.txt',fpl(i,:), 'Delimiter',' ', '-append', 'precision', '%.40f')
            else
                disp([num2str(i),': Good']);
            end
        else
            disp(['Drela: ',num2str(double(nint))]);
            disp(['Me: ',num2str(q_ind)]);
            disp(['Diff: ',num2str(double(d))]);
        end
    catch
        disp([num2str(i),': Error. FPL = [',num2str(fpl(i,1)),', ', num2str(fpl(i,2)), ', ', num2str(fpl(i,3)), ']'])
        dlmwrite('BadFPs.txt',fpl(i,:), 'Delimiter',' ', '-append', 'precision', '%.40f')
    end
end

%% Plot

% [hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCONTROL, matPLEX, [], vecUINF, matROTANG, [3 1 4 4], 'opengl');
% 
% % granularity = 0.0025;
% % x = [-0.05:granularity:0.05];
% % y = [0.01];
% % z = [-0.05:granularity:0.05];
% % 
% % granularity = 0.0025;
% % x = [0.225:granularity:0.275];
% % y = [0.25];
% % z = [-0.05:granularity:0.05];
% 
% % granularity = 0.05;
% % x = [-5:granularity:7];
% % y = [-5, matCENTER(:,2), 5];
% % z = [-5:granularity:5];
% 
% [X,Y,Z] = meshgrid(x,y,z);
% fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');
% [q_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCONTROL);
% figure(1);
% % q_ind = permute(q_ind,[3 2 1]);
% hold on
% quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3),'b')
% hold off







