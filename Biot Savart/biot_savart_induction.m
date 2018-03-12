clc
clear

matPOINTS(:,:,1) = [0 0 0];
matPOINTS(:,:,2) = [0 2 0];
matPOINTS(:,:,3) = [1 0.5 0];

matCOEFF = [1 1 1 1 1 1];

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

granularity = 0.20;
y = -0.2:granularity:2;
z = -.4:granularity:0.4;
x = -.4:granularity:1.3;

% granularity = 0.5;
% y = 1.8;
% z = -1:0.1:1;
% x = 2;

% granularity = 0.4;
% y = -9:granularity:9;
% z = -2:granularity:2;
% x = 4;

[X,Y,Z] = meshgrid(x,y,z);
fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');

len = size(fpg,1);
dvenum = ones(len, 1);
infl_glob = induc(repmat(matPOINTS,len,1,1), fpg);
q_ind = permute(sum(infl_glob.*repmat(reshape(matCOEFF(dvenum,:)',1,6,[]),3,1,1),2),[2 1 3]);
q_ind = reshape(permute(q_ind,[3 1 2]),[],3,1)./(-4*pi);

% tic
% for i = 1:size(fpg,1)
%   q_ind(i,:) = induc(matPOINTS, fpg(i,:), matCOEFF);
% end
% toc

fpg(any(isinf(q_ind),2),:) = [];
q_ind(any(isinf(q_ind),2),:) = [];

hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3))
hold off

function [infl_glob] = induc(matPOINTS, fpg)
xi_1 = matPOINTS(:,1,1);
xi_2 = matPOINTS(:,1,2);
xi_3 = matPOINTS(:,1,3);

eta_1 = matPOINTS(:,2,1);
eta_2 = matPOINTS(:,2,2);
eta_3 = matPOINTS(:,2,3);

xi_p = fpg(:,1);
eta_p = fpg(:,2);
zeta_p = fpg(:,3);

% A1 = matCOEFF(:,1);
% A2 = matCOEFF(:,2);
% B1 = matCOEFF(:,3);
% B2 = matCOEFF(:,4);
% C2 = matCOEFF(:,5);

syms x y
syms xi eta

% t2 = abs(zeta_p);
% t6 = y-eta_p;
% t3 = abs(t6);
% t8 = x-xi_p;
% t4 = abs(t8);
% t5 = t2.^2;
% t7 = t3.^2;
% t9 = t4.^2;
% t10 = t5+t7+t9;
% t11 = 1.0./t10.^(3.0./2.0);
% t12 = A1.*y;
% t13 = C2.*x;
% t14 = A2+t12+t13;
% t15 = C2.*y;
% t16 = B1.*x;
% t17 = B2+t15+t16;
% tr = [-t11.*t17.*zeta_p,-t11.*t14.*zeta_p,-t11.*(t6.*t14+t8.*t17)];

d_inf = [...
    0,...
    -eta.*zeta_p.*1.0./(abs(zeta_p).^2+abs(eta-eta_p).^2+abs(xi-xi_p).^2).^(3.0./2.0),...
    -eta.*(eta-eta_p).*1.0./(abs(zeta_p).^2+abs(eta-eta_p).^2+abs(xi-xi_p).^2).^(3.0./2.0), ...
    ...
    0,...
    -zeta_p.*1.0./(abs(zeta_p).^2+abs(eta-eta_p).^2+abs(xi-xi_p).^2).^(3.0./2.0), ...
    -(eta-eta_p).*1.0./(abs(zeta_p).^2+abs(eta-eta_p).^2+abs(xi-xi_p).^2).^(3.0./2.0), ...
    ...
    -xi.*zeta_p.*1.0./(abs(zeta_p).^2+abs(eta-eta_p).^2+abs(xi-xi_p).^2).^(3.0./2.0), ...
    0, ...
    -xi.*(xi-xi_p).*1.0./(abs(zeta_p).^2+abs(eta-eta_p).^2+abs(xi-xi_p).^2).^(3.0./2.0),...
    -zeta_p.*1.0./(abs(zeta_p).^2+abs(eta-eta_p).^2+abs(xi-xi_p).^2).^(3.0./2.0),...
    0.0,...
    -(xi-xi_p).*1.0./(abs(zeta_p).^2+abs(eta-eta_p).^2+abs(xi-xi_p).^2).^(3.0./2.0),...
    -eta.*zeta_p.*1.0./(abs(zeta_p).^2+abs(eta-eta_p).^2+abs(xi-xi_p).^2).^(3.0./2.0),-xi.*zeta_p.*1.0./(abs(zeta_p).^2+abs(eta-eta_p).^2+abs(xi-xi_p).^2).^(3.0./2.0),...
    -(xi.*(eta-eta_p)+eta.*(xi-xi_p)).*1.0./(abs(zeta_p).^2+abs(eta-eta_p).^2+abs(xi-xi_p).^2).^(3.0./2.0),...
    ...
    0.0,0.0,0.0];

d_inf = subs(d_inf, 'eta', 'y');
d_inf = subs(d_inf, 'xi', 'x');

eta_le = eta_2 + (x - xi_2).*((eta_3 - eta_2)./(xi_3 - xi_2));
eta_te = eta_1 + (x - xi_1).*((eta_3 - eta_1)./(xi_3 - xi_1));

count = 1;
h = waitbar(count/length(d_inf(:)),'Please wait...');
infl_glob = zeros(3,6,size(fpg,1));
for i = 1:size(d_inf,1)
    for j = 1:size(d_inf,2)
        if isempty(find(d_inf(i,j)==0,1))
            infl_glob(mod(j-1,3) + 1, floor(j/3) + 1, i) = integral2(matlabFunction(d_inf(i,j)), xi_1(i), xi_3(i), matlabFunction(eta_te(i)), matlabFunction(eta_le(i)));
        end
        waitbar(count/length(d_inf(:)), h)
        count = count + 1;
    end
end

close(h);

end
