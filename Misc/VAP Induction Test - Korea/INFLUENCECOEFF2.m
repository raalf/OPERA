function [D] = INFLUENCECOEFF2(testDVE, DVECT, FP)

num = 1;

numsides = size(testDVE,1); %number of sides
Xi_1 = testDVE( [1:1:numsides],:); % Vertex 1
Xi_2 = testDVE( [[2:1:numsides] , 1],:); % Vertex 2

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

%k, normal to triangle formed with FP
K = cross(R_i,R_j);
ki = sqrt((K(:,1).^2 + K(:,2).^2 + K(:,3).^2));
K = K./repmat(ki,[1,size(s_i,2)]) ;
[idx6,~] = find(~ki);
K(idx6,:) = zeros(length(idx6),3);

N_2D = cross(K,S_i);

%field point in local frame
x_prime = dot(-R_i,S_i,2);

y_prime = dot(-R_i,N_2D,2);

r_i = sqrt(x_prime.^2 + y_prime.^2);
r_j = sqrt((l_i-x_prime).^2 + y_prime.^2);

%unknown variables
e_n = e_zeta;
temp = dot(repmat(e_zeta,[size(Xi_1,1),1]),R_i);
e_n(temp<0) = -e_zeta(temp<0);

e = dot(repmat(e_n,[size(Xi_1,1),1]), R_i,2);
b_i = dot(cross(repmat(e_zeta,[size(Xi_1,1),1]),R_i,2),S_i,2);

%edge in local frame
si_xsi = dot(S_i,repmat(e_xsi,[ size(s_i,1),1]),2); 
si_eta = dot(S_i,repmat(e_eta,[ size(s_i,1),1]),2);

%angle check
F = (( y_prime.^2 + e.*r_i) ./ (r_i + e)) .^2 + ((y_prime.^2 + e.*r_j)./(r_j + e)).^2 - y_prime.^2;

% WHAT DO WE DO IT Y_PRIME^2 < e^2 ??????????????????????????????????????????????????????????????????????

H = (sqrt(y_prime.^2 - e.^2) .* (y_prime.^2 .*l_i + e.*(l_i - x_prime).*r_i + e.*x_prime.*r_j ))./...
    (y_prime.^2.*(r_i + e).*(r_j+e));

[indx1,n] = find(F>0);
[indx2,n] = find(F<=0);

beta(indx1,1) = asin(H(indx1))';
beta(indx2,1) = pi-asin(H(indx2))';
if min(F) <= 0
   disp('Beta conditional') 
end

%unknown
E = log10((r_j + l_i - x_prime)./ (r_i - x_prime) );
% r_j
% l_i
% x_prime
% r_i

%solve K1&K2
% if abs(y_prime) - e <= 1e-10
test_2 = y_prime - e;
[indx3,~] = find(test_2 == 0);

K1(indx3,1) = E(indx3) - (x_prime(indx3)./(r_i(indx3) + e(indx3))) - ((l_i(indx3) - x_prime(indx3))./(r_j(indx3) + e(indx3)));
% WHAT TO DO IF Y_PRIME < e??????
% K1(test_2 ~= 0,1) = E(test_2 ~= 0) - (e(test_2 ~= 0)./sqrt(y_prime(test_2 ~= 0).^2 - e(test_2 ~= 0))).^2.*beta(test_2 ~= 0);
K1(test_2 ~= 0,1) = E(test_2 ~= 0) - (e(test_2 ~= 0)./sqrt(y_prime(test_2 ~= 0).^2 - e(test_2 ~= 0).^2)).*beta(test_2 ~= 0);
% K1
% E
% x_prime
% y_prime
% e
% beta

K2 = (1/2)*((l_i - x_prime).*r_j + x_prime.*r_i + y_prime.^2.*E);

K1(b_i == 0) = 0;
K1(isinf(E)) = 0;
K2(y_prime == 0) = (1/2)*((l_i(y_prime == 0) - x_prime(y_prime == 0)).*r_j(y_prime == 0) + x_prime(y_prime == 0).*r_i(y_prime == 0)); 
K2(isinf(E)) = (1/2)*((l_i(isinf(E)) - x_prime(isinf(E))).*r_j(isinf(E)) + x_prime(isinf(E)).*r_i(isinf(E)));

%solve M,N,P 
M = sum(b_i.*K1);
N = sum(si_xsi.*K2);
P = sum(si_eta.*K2);

% M = sum(b_i.*K1);
% N = sum(si_eta.*K2);
% P = sum(si_xsi.*K2);

% test_var = cross(e_zeta,R_i(1,:));
% test_var = cross(K(1,:),R_i(1,:))
% test_var2 = cross(R_i(1,:),S_i(1,:));

% hFig2 = figure(2);
% clf(2)
% patch(testDVE(:,1),testDVE(:,2),testDVE(:,3),'r','LineWidth',2)
% alpha(0.5);
% hold on
% scatter3(FP(1),FP(2),FP(3),100,'ok','filled');
% %r_i
% plot3([Xi_1(1,1) Xi_1(1,1)-R_i(1,1)],[Xi_1(1,2) Xi_1(1,2)-R_i(1,2)],[Xi_1(1,3) Xi_1(1,3)-R_i(1,3)],'k','LineWidth',2)
% %r_j
% plot3([Xi_2(1,1) Xi_2(1,1)-R_j(1,1)],[Xi_2(1,2) Xi_2(1,2)-R_j(1,2)],[Xi_2(1,3) Xi_2(1,3)-R_j(1,3)],'b','LineWidth',2)
% %k
% % plot3([Xi_1(1,1) Xi_1(1,1)+K(1,1)],[Xi_1(1,2) Xi_1(1,2)+K(1,2)],[Xi_1(1,3) Xi_1(1,3)+K(1,3)],'-m','LineWidth',2)
% %s_i
% % plot3([Xi_1(1,1) Xi_1(1,1)+S_i(1,1)],[Xi_1(1,2) Xi_1(1,2)+S_i(1,2)],[Xi_1(1,3) Xi_1(1,3)+S_i(1,3)],'-g','LineWidth',2)
% 
% % plot3([Xi_1(1,1) Xi_1(1,1)+test_var(1)],[Xi_1(1,2) Xi_1(1,2)+test_var(2)],[Xi_1(1,3) Xi_1(1,3)+test_var(3)],'-m','LineWidth',2)
% 
% % plot3([Xi_1(1,1) Xi_1(1,1)+test_var2(1)],[Xi_1(1,2) Xi_1(1,2)+test_var2(2)],[Xi_1(1,3) Xi_1(1,3)+test_var2(3)],'--b','LineWidth',2)
% 
% hold off
% grid on
% axis equal
% box on

D(1,:) = [-P*e_eta(1) M*e_zeta(1) N*e_xsi(1) M*e_zeta(1) 0];
D(2,:) = [-P*e_eta(2) M*e_zeta(2) N*e_xsi(2) M*e_zeta(2) 0];
D(3,:) = [-P*e_eta(3) M*e_zeta(3) N*e_xsi(3) M*e_zeta(3) 0];

% D(1,:) = [P*e_xsi(1) M*e_zeta(1) -N*e_eta(1) M*e_zeta(1) 0];
% D(2,:) = [P*e_xsi(2) M*e_zeta(2) -N*e_eta(2) M*e_zeta(2) 0];
% D(3,:) = [P*e_xsi(3) M*e_zeta(3) -N*e_eta(3) M*e_zeta(3) 0];

% D(1,:) = [-P 0 0 0 0];
% D(2,:) = [0 0 N 0 0];
% D(3,:) = [0 M 0 M 0];

end

