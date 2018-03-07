clc
clear

matPOINTS(:,:,1) = [0 0 0];
matPOINTS(:,:,2) = [0 2 0];
matPOINTS(:,:,3) = [1 0.5 0];

matCOEFF = [0 1 0 0 0 0];

hFig203 = figure(203);
clf(203);
patch(reshape(matPOINTS(:,1,:),[],1,1),reshape(matPOINTS(:,2,:),[],1,1),reshape(matPOINTS(:,3,:),[],1,1),'FaceColor','red','FaceAlpha',0.3,'LineWidth',2)
axis equal
axis tight
grid minor
box on
xlabel('\xi-Direction','FontSize',15);
ylabel('\eta-Direction','FontSize',15);
zlabel('\zeta-Direction','FontSize',15);

% [TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
%     matELOC, matPLEX, matDVECT, matALIGN, matVATT, matVNORM, matCENTER, matROTANG] = fcnTRIANG(matPOINTS);

% fpg = [2 5 1];

% granularity = 0.1;
% y = -0.2:granularity:2;
% z = .4:granularity:0.4;
% x = -.4:granularity:1.3;

% granularity = 0.5;
% y = 1.8;
% z = -1:0.1:1;
% x = 2;

granularity = 0.4;
y = -9:granularity:9;
z = -2:granularity:2;
x = 4;

[X,Y,Z] = meshgrid(x,y,z);
fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');

tic
for i = 1:size(fpg,1)
  q_ind(i,:) = induc(matPOINTS, fpg(i,:), matCOEFF);
end
toc

fpg(any(isinf(q_ind),2),:) = [];
q_ind(any(isinf(q_ind),2),:) = [];

hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3))
hold off

function [q_ind] = induc(matPOINTS, fpg, matCOEFF)
xi_1 = matPOINTS(1,1,1);
xi_2 = matPOINTS(1,1,2);
xi_3 = matPOINTS(1,1,3);

eta_1 = matPOINTS(1,2,1);
eta_2 = matPOINTS(1,2,2);
eta_3 = matPOINTS(1,2,3);

xi_p = fpg(1);
eta_p = fpg(2);
zeta_p = fpg(3);

A1 = matCOEFF(1);
A2 = matCOEFF(2);
B1 = matCOEFF(3);
B2 = matCOEFF(4);
C2 = matCOEFF(5);

q_xi = w_xi(B1,B2,C2,eta_1,eta_2,eta_3,xi_1,xi_2,xi_3,xi_p,zeta_p);
q_zeta = w_zeta(A1,A2,B1,B2,C2,eta_1,eta_2,eta_3,eta_p,xi_1,xi_3,xi_p,zeta_p);

q_ind = [q_xi q_xi.*0 q_zeta];

end
