clc
clear

matPOINTS(:,:,1) = [0 0 0];
matPOINTS(:,:,2) = [0 2 0];
matPOINTS(:,:,3) = [1 0.5 0];

matCOEFF = [0 0 1 1 0 0];

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

granularity = 0.20;
y = -0.2:granularity:2;
z = -.4:granularity:0.4;
x = -.4:granularity:1.3;

% granularity = 0.025;
% y = 0.6;
% z = -0.2:granularity:0.2;
% x = 0.4;

% % granularity = 0.5;
% % y = 1.8;
% % z = -1:0.1:1;
% % x = 2;

% granularity = 0.4;
% y = -9:granularity:9;
% z = -2:granularity:2;
% x = 4;
 
[X,Y,Z] = meshgrid(x,y,z);
fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');

% fpg = [0.4 1.4 0];

len = size(fpg,1);
dvenum = ones(len, 1);
% infl_glob = induc(repmat(matPOINTS,len,1,1), fpg);
load('matlab.mat','infl_glob');

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

function [infl_loc] = induc(matPOINTS, fpl)

xi_1 = matPOINTS(:,1,1);
xi_2 = matPOINTS(:,1,2);
xi_3 = matPOINTS(:,1,3);

eta_1 = matPOINTS(:,2,1);
eta_2 = matPOINTS(:,2,2);
eta_3 = matPOINTS(:,2,3);

% Checking which elements are on the element, moving them off by a small
% amount
margin = 1e-5;
le_eta = eta_2 + (fpl(:,1) - xi_2).*((eta_3 - eta_2)./(xi_3 - xi_2));
te_eta = eta_1 + (fpl(:,1) - xi_1).*((eta_3 - eta_1)./(xi_3 - xi_1));
idx_on_element = fpl(:,2) >= te_eta & fpl(:,2) <= le_eta & fpl(:,1) >= xi_1 & fpl(:,1) <= xi_3 & abs(fpl(:,3)) <= margin;
idx_on_edge = idx_on_element & ((le_eta - margin <= fpl(:,2) & fpl(:,2) <= le_eta + margin) | (te_eta - margin <= fpl(:,2) & fpl(:,2) <= te_eta + margin));

fpl(idx_on_element,3) = fpl(idx_on_element,3) + 100*margin;
fpl(idx_on_edge,2) = fpl(idx_on_edge,2) + 100*margin;

xi_p = fpl(:,1);
eta_p = fpl(:,2);
zeta_p = fpl(:,3);

syms x y
syms xi eta

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
infl_loc = zeros(3,6,size(fpl,1));
for i = 1:size(d_inf,1)
    for j = 1:size(d_inf,2)
        if isempty(find(d_inf(i,j)==0,1))
            infl_loc(mod(j-1,3) + 1, floor(j/3) + 1, i) = integral2(matlabFunction(d_inf(i,j)), xi_1(i), xi_3(i), matlabFunction(eta_te(i)), matlabFunction(eta_le(i)));
        end
        waitbar(count/length(d_inf(:)), h)
        count = count + 1;
    end
end
close(h);

% Only using the tangential velocities for points on the element
infl_loc(1:2,:,idx_on_element) = infl_loc(1:2,:,idx_on_element).*0;

end
