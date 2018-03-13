clc
clear

%% Defining element
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

%% Points to influence on 
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

% fpg = [0.4 0.6 0.2];

%% Calculating induced velocities
len = size(fpg,1);
dvenum = ones(len, 1);
infl_glob = induc(repmat(matPOINTS,len,1,1), fpg);
% load('matlab.mat')

q_ind = permute(sum(infl_glob.*repmat(reshape(matCOEFF(dvenum,:)',1,6,[]),3,1,1),2),[2 1 3]);
q_ind = reshape(permute(q_ind,[3 1 2]),[],3,1)./(-4*pi);

hold on
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3))
hold off

%% induction function
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
idx_on_element = fpl(:,2) >= te_eta - margin & fpl(:,2) <= le_eta + margin & fpl(:,1) >= xi_1 & fpl(:,1) <= xi_3 & abs(fpl(:,3)) <= margin;
idx_on_edge = idx_on_element & ((xi_1 - margin <= fpl(:,1) & fpl(:,1) <= xi_1 + margin) | (le_eta - margin <= fpl(:,2) & fpl(:,2) <= le_eta + margin) | (te_eta - margin <= fpl(:,2) & fpl(:,2) <= te_eta + margin));

fpl(idx_on_element,3) = fpl(idx_on_element,3) + 10*margin;
fpl(idx_on_edge,2) = fpl(idx_on_edge,2) + 1000*margin;

xi_p = fpl(:,1);
eta_p = fpl(:,2);
zeta_p = fpl(:,3);

len = size(fpl,1);

count = 1;
h = waitbar(count/len,'Please wait...');
infl_loc = zeros(3,6,len);

for i = 1:len
    
    le_eta = @(x) eta_2(i) + (x - xi_2(i)).*((eta_3(i) - eta_2(i))./(xi_3(i) - xi_2(i)));
    te_eta = @(x) eta_1(i) + (x - xi_1(i)).*((eta_3(i) - eta_1(i))./(xi_3(i) - xi_1(i)));
    
    for j = 1:15
        % a1
        if j == 2
            term = @(x,y) -y.*zeta_p(i).*1.0./(abs(zeta_p(i)).^2+abs(y-eta_p(i)).^2+abs(x-xi_p(i)).^2).^(3.0./2.0);
            infl_loc(2,1,i) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta);
        elseif j == 3
            term = @(x,y) -y.*(y-eta_p(i)).*1.0./(abs(zeta_p(i)).^2+abs(y-eta_p(i)).^2+abs(x-xi_p(i)).^2).^(3.0./2.0);
            infl_loc(3,1,i) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta); 
        % a2    
        elseif j == 5
            term = @(x,y) -zeta_p(i).*1.0./(abs(zeta_p(i)).^2+abs(y-eta_p(i)).^2+abs(x-xi_p(i)).^2).^(3.0./2.0);
            infl_loc(2,2,i) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta);               
        elseif j == 6
            term = @(x,y) -(y-eta_p(i)).*1.0./(abs(zeta_p(i)).^2+abs(y-eta_p(i)).^2+abs(x-xi_p(i)).^2).^(3.0./2.0);
            infl_loc(3,2,i) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta);
        % b1    
        elseif j == 7
            term = @(x,y) -x.*zeta_p(i).*1.0./(abs(zeta_p(i)).^2+abs(y-eta_p(i)).^2+abs(x-xi_p(i)).^2).^(3.0./2.0);
            infl_loc(1,3,i) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta);             
        elseif j == 9
            term = @(x,y) -x.*(x-xi_p(i)).*1.0./(abs(zeta_p(i)).^2+abs(y-eta_p(i)).^2+abs(x-xi_p(i)).^2).^(3.0./2.0);
            infl_loc(3,3,i) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta); 
        % b2    
        elseif j == 10
            term = @(x,y) -zeta_p(i).*1.0./(abs(zeta_p(i)).^2+abs(y-eta_p(i)).^2+abs(x-xi_p(i)).^2).^(3.0./2.0);
            infl_loc(1,4,i) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta);             
        elseif j == 12
            term = @(x,y) -(x-xi_p(i)).*1.0./(abs(zeta_p(i)).^2+abs(y-eta_p(i)).^2+abs(x-xi_p(i)).^2).^(3.0./2.0);
            infl_loc(3,4,i) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta);  
        % c2    
        elseif j == 13
            term = @(x,y) -y.*zeta_p(i).*1.0./(abs(zeta_p(i)).^2+abs(y-eta_p(i)).^2+abs(x-xi_p(i)).^2).^(3.0./2.0);
            infl_loc(1,5,i) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta);             
        elseif j == 14
            term = @(x,y) -x.*zeta_p(i).*1.0./(abs(zeta_p(i)).^2+abs(y-eta_p(i)).^2+abs(x-xi_p(i)).^2).^(3.0./2.0);
            infl_loc(2,5,i) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta);   
        elseif j == 15
            term = @(x,y) -(x.*(y-eta_p(i))+y.*(x-xi_p(i))).*1.0./(abs(zeta_p(i)).^2+abs(y-eta_p(i)).^2+abs(x-xi_p(i)).^2).^(3.0./2.0);
            infl_loc(3,5,i) = integral2(term, xi_1(i), xi_3(i), te_eta, le_eta);
        end        
    end
    
    waitbar(count/len, h)
    count = count + 1;
    
end
close(h);

% Only using the tangential velocities for points on the element
infl_loc(1:2,:,idx_on_element) = infl_loc(1:2,:,idx_on_element).*0;
% Zero it all for singular edges
infl_loc(:,:,idx_on_edge) = infl_loc(:,:,idx_on_edge).*0;

end
