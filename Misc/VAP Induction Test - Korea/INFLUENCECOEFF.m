function [M,N,P]=INFLUENCECOEFF(testDVE,DNORM,DVECT,FP)
num = 1;

numsides = size(testDVE,1); %number of sides
Xi_1 = testDVE( [1:1:numsides],:); % Vertex 1
Xi_2 = testDVE( [[2:1:numsides] , 1],:); % Vertex 2

%% 

%element local coords in global
e_eta = DVECT(num,:,1);
e_xsi = DVECT(num,:,2);
e_zeta = DVECT(num,:,3); %should be element normal

%vectors from vertex to field point
R_i = Xi_1 - repmat(FP,[size(Xi_1,1),1]);
R_j = Xi_2 - repmat(FP,[size(Xi_1,1),1]); % R_(i+1)

%vector of edge
s_i = Xi_2 - Xi_1;
l_i = sqrt((s_i(:,1).^2 + s_i(:,2).^2 + s_i(:,3).^2));
S_i = s_i./repmat(l_i,[1,size(s_i,2)]); 

%% THEIR WAY
% S_i = edge of dve
% e_zeta = element normal
% K = S_i x Normal
% K = cross(R_i,R_j);

K = cross(S_i,repmat(e_zeta,[size(Xi_1,1),1]),2);

% p_norm = cross(R_i,R_j);
% p_norm_norm = sqrt(sum(abs(p_norm).^2,2));
% K = K - (dot(K,p_norm,2)./(p_norm_norm).^2).*p_norm;
% K = K./(sqrt(sum(abs(K).^2,2)));

%field point in local frame
x_prime = dot(-R_i,S_i,2);
% y_prime = dot(cross(R_i,S_i,2),repmat(e_zeta,[size(Xi_1,1),1]),2);%!!!!!!
y_prime = dot(cross(R_i,S_i,2),K,2);

r_i = sqrt(x_prime.^2 + y_prime.^2);
r_j = sqrt((l_i-x_prime).^2 + y_prime.^2);

%unknown variables
e = abs(dot(repmat(e_zeta,[size(Xi_1,1),1]), R_i,2));
b_i = dot(cross(repmat(e_zeta,[size(Xi_1,1),1]),R_i,2),S_i,2);

%edge in local frame
si_xsi = dot(S_i,repmat(e_xsi,[ size(s_i,1),1]),2); 
si_eta = dot(S_i,repmat(e_eta,[ size(s_i,1),1]),2);

%angle check
F = (( y_prime.^2 + e.*r_i) ./ (r_i + e)) .^2 + ((y_prime.^2 + e.*r_j)./(r_j + e)).^2 - y_prime.^2;
H = (sqrt(y_prime.^2 - e.^2) .* (y_prime.^2 .*l_i + e.*(l_i - x_prime).*r_i + e.*x_prime.*r_j ))./...
    (y_prime.^2.*(r_i + e).*(r_j+e));
[indx1,n] = find(F>0);
[indx2,n] = find(F<=0);

beta(indx1,:) = asin(H(indx1));
beta(indx2) = pi-asin(H(indx2));
clear n indx1 indx2

%unknown
E = log10((r_j + l_i - x_prime)./ (r_i - x_prime) );

%solve K1&K2
K1 = E- (e./sqrt(y_prime.^2 - e.^2)).*beta
K2 = (1/2)*((l_i - x_prime).*r_j + x_prime.*r_i + y_prime.^2.*E)

%solve M,N,P 
M = b_i(1)*K1(1) + b_i(2)*K1(2) + b_i(3)*K1(3);
N = si_xsi(1)*K2(1) + si_xsi(2)*K2(2) + si_xsi(3)*K2(3);
P = si_eta(1)*K2(1) + si_eta(2)*K2(2) + si_eta(3)*K2(3);

test_var = cross(e_zeta,R_i(1,:));
% test_var = cross(K(1,:),R_i(1,:))
test_var2 = cross(R_i(1,:),S_i(1,:));

hFig1 = figure(1);
clf(1)
patch(testDVE(:,1),testDVE(:,2),testDVE(:,3),'r','LineWidth',2)
alpha(0.5);
hold on
scatter3(FP(1),FP(2),FP(3),100,'ok','filled');
%r_i
plot3([Xi_1(1,1) Xi_1(1,1)-R_i(1,1)],[Xi_1(1,2) Xi_1(1,2)-R_i(1,2)],[Xi_1(1,3) Xi_1(1,3)-R_i(1,3)],'k','LineWidth',2)
%r_j
plot3([Xi_2(1,1) Xi_2(1,1)-R_j(1,1)],[Xi_2(1,2) Xi_2(1,2)-R_j(1,2)],[Xi_2(1,3) Xi_2(1,3)-R_j(1,3)],'b','LineWidth',2)
%k
plot3([Xi_1(1,1) Xi_1(1,1)+K(1,1)],[Xi_1(1,2) Xi_1(1,2)+K(1,2)],[Xi_1(1,3) Xi_1(1,3)+K(1,3)],'-m','LineWidth',2)
%s_i
plot3([Xi_1(1,1) Xi_1(1,1)+S_i(1,1)],[Xi_1(1,2) Xi_1(1,2)+S_i(1,2)],[Xi_1(1,3) Xi_1(1,3)+S_i(1,3)],'-g','LineWidth',2)

plot3([Xi_1(1,1) Xi_1(1,1)+test_var(1)],[Xi_1(1,2) Xi_1(1,2)+test_var(2)],[Xi_1(1,3) Xi_1(1,3)+test_var(3)],'-y','LineWidth',2)

plot3([Xi_1(1,1) Xi_1(1,1)+test_var2(1)],[Xi_1(1,2) Xi_1(1,2)+test_var2(2)],[Xi_1(1,3) Xi_1(1,3)+test_var2(3)],'--b','LineWidth',2)

hold off
grid on
axis equal
box on


%plotting
% figure(1)
% clf(1)
% hold on
% patch(testDVE(:,1),testDVE(:,2),testDVE(:,3),'b')
% plot3(FP(1),FP(2),FP(3),'*k')

% MY WAY (not projecting the distance vectors) - this is likely not right
% x_prime = PLEX(2,1) % Local eta coord of FP
% y_prime = PLEX(2,2) % Local xi coord of FP
% 
% l_i = PLEX(3,1) % Local eta coord of Vertex 2 point
% 
% r_i = sqrt(x_prime^2 + y_prime^2)
% r_j = sqrt((l_i-x_prime)^2 + y_prime^2)

% theta_i = atan((abs(y_prime)*l_i)/(r_i^2 - l_i*x_prime));
% % theta_i = atan2((abs(y_prime)*l_i),(r_i^2 - l_i*x_prime))
% I_1 = ((l_i - x_prime)*log10(r_j^2)) + (x_prime*log10(r_i^2)) - (2*l_i) + (2*abs(y_prime)*theta_i);
% I_2 = 0.5*(r_j^2 * log10(r_j^2) - r_i^2*log10(r_i^2)) - 0.5*l_i^2 + l_i*x_prime + x_prime*I_1;
% 
% syms omega
% syms grd_omega
% I = 0.5.*N_i.*(omega(l_i + I_1) + (dot(grd_omega,S_i)*(0.5*l_i^2 + I_2))) - 0.25*grd_omega(dot(N_i,R_i))*I_1
