clear
% clc

%% Preamble
set(groot,'defaultFigureCreateFcn','addToolbarExplorationButtons(gcf)')
set(groot,'defaultAxesCreateFcn','set(get(gca,''Toolbar''),''Visible'',''off'')')

strFILE = 'inputs/2dve.dat';

% [~, strATYPE, vecSYM, ~, ~, ~, valALPHA, valBETA, ~, ~] = fcnOPREAD(strFILE);
strATYPE = 'what';
valALPHA = 0;
valBETA = 0;
vecSYM = [];


% matPOINTS(:,:,1) = [0 0 0];
% matPOINTS(:,:,2) = [0  0.5 0];
% matPOINTS(:,:,3) = [0.5  0.5 0];

matPOINTS(:,:,1) = [1 0 0];
matPOINTS(:,:,2) = [0  0 0];
matPOINTS(:,:,3) = [0.3  0.5 0];


[TR, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
    matELOC, matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG, matCONTROL] = fcnTRIANG(matPOINTS, 'SURFACE', []);

vecUINF = fcnUINFWING(valALPHA, 0);

matCOEFF = [0 0 -2 0.17 0 0.056];

%% fpl
% fpg = [-0 0.2 0]

% % Wide
% granularity = 0.1;
% x = [-5:granularity:7];
% y = [-5:granularity:5];
% z = [-5:granularity:5];
% [X,Y,Z] = meshgrid(x,y,z);
% fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');

% % Narrow
% granularity = 0.05;
% x = [-0.05:granularity:0.55];
% y = [-0.05:granularity:0.55];
% z = [-0.1:granularity:0.1];
% [X,Y,Z] = meshgrid(x,y,z);
% fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');

% fpg = [0.5 0.0000000000000000000000000000000000000000 0.001]

% fpl = fcnGLOBSTAR(fpg - matCENTER, matROTANG);
% len = size(fpl,1);

%% Coefficients
% syms xi eta real
% mu = 0.5.*matCOEFF(1).*eta.^2 + matCOEFF(2).*eta + 0.5.*matCOEFF(3).*xi.^2 + matCOEFF(4).*xi + matCOEFF(5);
% 
% s = [xi eta xi.*0];
% n = matDVECT(:,:,3);
% 
% xi1 = matPLEX(1,1);
% xi2 = matPLEX(2,1);
% xi3 = matPLEX(3,1);
% eta1 = matPLEX(1,2);
% eta2 = matPLEX(2,2);
% eta3 = matPLEX(3,2);
% 
% le = (((eta3 - eta2).*xi)./(xi3 - xi2)) + eta2 - ((xi2.*(eta3 - eta2))./(xi3 - xi2));
% te = (((eta3 - eta1).*xi)./(xi3 - xi1)) + eta1 - ((xi1.*(eta3 - eta1))./(xi3 - xi1));
% 
% infl_loc = fcnHDVEIND_DB(ones(len,1), ones(len,1), fpg, matPLEX, matROTANG, matCENTER);
% q_ind = permute(sum(infl_loc.*repmat(reshape(matCOEFF(1,:)',1,5,[]),3,1,1),2),[2 1 3]);
% 
% if len > 1
%     delete('BadFPs.txt');
% end
% 
% % for i = 1:size(fpl,1)
% parfor i = 1:size(fpl,1)
%    
%     infl_loc = fcnHDVEIND_DB(1, 1, fpg(i,:), matPLEX, matROTANG, matCENTER);
%     q_ind = permute(sum(infl_loc.*repmat(reshape(matCOEFF(1,:)',1,5,[]),3,1,1),2),[2 1 3]);
%      
%     try
%         r = fpl(i,:);
%         term = repmat(mu,1,3).*((n./(sqrt(sum((r - s).^2,2)).^3)) - dot(3.*repmat(n,length(xi(:)),1), r - s, 2).*((r - s)./(sqrt(sum((r - s).^2,2)).^5)));
%         nint = vpaintegral(int(term, eta, le, te), xi, xi1, xi3);
%         
%         d = abs((nint - q_ind)./q_ind).*100;
%         
%         if len > 1
%             if max(real(d)) > 2e-0 && ~any(isnan(d))
%                 disp([num2str(i),': Bad. FPG = [',num2str(fpg(i,1)),', ', num2str(fpg(i,2)), ', ', num2str(fpg(i,3)), ']'])
%                 dlmwrite('BadFPs.txt',fpg(i,:), 'Delimiter',' ', '-append', 'precision', '%.40f')
%             else
%                 disp([num2str(i),': Good']);
%             end
%         else
%             disp(['Drela: ',num2str(double(nint))]);
%             disp(['Me: ',num2str(q_ind)]);
%             disp(['Diff: ',num2str(double(d))]);
%         end
%     catch
%         disp([num2str(i),': VPA Failed. FPG = [',num2str(fpg(i,1)),', ', num2str(fpg(i,2)), ', ', num2str(fpg(i,3)), ']'])
%         %         dlmwrite('BadFPs.txt',fpg(i,:), 'Delimiter',' ', '-append', 'precision', '%.40f')
%     end
% end

%% Plot

[hFig1] = fcnPLOTBODY(1, matDVE, valNELE, matVLST, matELST, matDVECT, matCONTROL, matPLEX, [], vecUINF, matROTANG, [3 1 4 4], 'opengl');
fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, [], matROTANG, 'r', 10);

% granularity = 0.000125;
% y = [-0.025:granularity:0.025];
% x = [0.2];
% z = [-0.025:granularity:0.025].*0;

% granularity = 0.00125;
granularity = 0.1
x = [0:granularity:0.2];
y = [0.2];
z = [-0.05 0 0.05];
% z = 0;

% granularity = 0.00125;
% % granularity = 0.1
% x = [0:granularity:0.2];
% y = [0.2];
% % z = [-0.05 0 0.05];
% z = 0;

% granularity = 0.0001;
% x = [0.225:granularity:0.275];
% y = [0.25];
% z = [-0.05:granularity:0.05].*0;

% granularity = 0.05;
% x = [-5:granularity:7];
% y = [-5, matCENTER(:,2), 5];
% z = [-5:granularity:5];

% granularity = 0.05;
% x = 0.2;
% y = [-0.1:granularity:0.1];
% z = 0.0002;

[X,Y,Z] = meshgrid(x,y,z);
fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');
% fpg = [0.5 0 0; 0  0 0; 0  0.5 0];
% fpg = [0 0.2 0];
% fpg = [0.2 0 0];
[q_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCONTROL);


%%
len = size(fpg,1);

dve2 = 1;
p1 = fcnGLOBSTAR(matVLST(matELST(1,1),:) - matCENTER(dve2,:), matROTANG(dve2,:));
p2 = fcnGLOBSTAR(matVLST(matELST(1,2),:) - matCENTER(dve2,:), matROTANG(dve2,:));

midpoint = mean([p1; p2],1);
vec = p2 - p1;
hspan = abs(p1(1,1) - p2(1,1))./2;
vec = vec./sqrt(sum(vec.^2,2));
phi = acos(dot([1 0 0], vec, 2));
if phi > pi/2
    phi = phi - pi;
end

xi = [linspace(p1(1,1), p2(1,1), 100)]';
eta = [linspace(p1(1,2), p2(1,2), 100)]';       
circ = 0.5.*matCOEFF(dve2,1).*eta.^2 + matCOEFF(dve2,2).*eta + 0.5.*matCOEFF(dve2,3).*xi.^2 + matCOEFF(dve2,4).*xi + matCOEFF(dve2,5).*xi.*eta + matCOEFF(dve2,6);
% circ = 0.5.*matCOEFF(dve2,3).*xi.^2 + matCOEFF(dve2,4).*xi + matCOEFF(dve2,5).*xi.*eta + matCOEFF(dve2,6);

if p2(1,1) > p1(1,1)
    tmp6 = [linspace(-hspan, hspan, 100)]';
else
    tmp6 = [linspace(hspan, -hspan, 100)]';
end     
fp_0 = fcnGLOBSTAR(fpg - matCENTER(dve2,:), matROTANG(dve2,:)) - midpoint;
[aloc, bloc, cloc] = fcnBOUNDIND(repmat(hspan,len,1), repmat(phi,len,1), [-fp_0(:,2) fp_0(:,1) fp_0(:,3)]);
%         aloc = [aloc(:,2), -aloc(:,1), aloc(:,3)];
%         bloc = [bloc(:,2), -bloc(:,1), bloc(:,3)];
%         cloc = [cloc(:,2), -cloc(:,1), cloc(:,3)]; 
        
D = [aloc bloc cloc]; % C is higher order
D = reshape(reshape(D', 1, 9, []), 3, 3, len);     
              
coeff = polyfit(tmp6, circ, 2);
coeff = fliplr(coeff);
w_ind = permute(sum(D.*repmat(reshape(coeff',1,3,[]),3,1,1),2),[2 1 3]);
w_ind = reshape(permute(w_ind,[3 1 2]),[],3,1)./(-4*pi);

% Out of filament local to element local
% w_ind = fcnSTARGLOB(w_ind, [phi.*0, phi.*0, phi]);

q_ind2 = q_ind + w_ind;
%%
figure(1);
% q_ind = permute(q_ind,[3 2 1]);
hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3),'b')
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), w_ind(:,1), w_ind(:,2), w_ind(:,3),'r')
% quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind2(:,1), q_ind2(:,2), q_ind2(:,3),'b')
hold off
setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera') ;


figure(10);
clf(10);
plot(fpg(:,1), q_ind(:,3), '-k');
hold on
plot(fpg(:,1), q_ind2(:,3), '--r');
hold off
grid minor
box on
axis tight
% set(gca, 'YScale', 'log')




