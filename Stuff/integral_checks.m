clc
clear
% warning off
tol = 1e-2;

len = 10;
%% 
c = linspace(-5,5,len)';
n = linspace(0,1,len)';

[C,N] = meshgrid(c,n);
fpg = unique([reshape(C,[],1) reshape(N,[],1)],'rows');
C = fpg(:,1);
N = fpg(:,2);
len = size(fpg,1);

alpha = zeros(len, 1) + 0.00001;
L = C;
u = N.^2 + alpha;
T = 2.*L.*N;
S = C.^2 + 1;

F_lim(:,:,1) = zeros(len, 1) + -1;
F_lim(:,:,2) = zeros(len, 1) + 1;

% F_lim(:,:,1) = zeros(len, 1) - 1;
% F_lim(:,:,2) = zeros(len, 1) - 0;

% F_lim(:,:,1) = zeros(len, 1) + 1;
% F_lim(:,:,2) = zeros(len, 1) + 2;

%%
F_lim = repmat(F_lim, 1, 1, 1);
S = repmat(S, 1, 2);
T = repmat(T, 1, 2);
u = repmat(u, 1, 2);

tol_alpha = 1e-6;
tol_S = 1e-3;
tol_T = 1e-3;
tol_u = 1e-2;
tol_ualpha = 1e-5;

idx = [];
F_lim = repmat(F_lim,1,2,1);
alpha = repmat(alpha,1,2,1);

% idx(:,:,1) = true(size(S));
idx(:,:,1:20) = false([size(S) 20]);

% idx(:,:,1)  = abs(S - 1) >  tol_S & abs(T) >  tol_T & abs(u) >  tol_u &  alpha >  tol_alpha;
% idx(:,:,2)  = abs(S - 1) >  tol_S & abs(T) >  tol_T & abs(u) >  tol_u &  alpha <= tol_alpha;
% idx(:,:,3)  = abs(S - 1) >  tol_S & abs(T) >  tol_T & abs(u) <= tol_u &  alpha >  tol_alpha;
% idx(:,:,4)  = abs(S - 1) >  tol_S & abs(T) >  tol_T & abs(u) <= tol_u &  alpha <= tol_alpha;
% idx(:,:,5)  = abs(S - 1) >  tol_S & abs(T) <= tol_T & abs(u) >  tol_u &  alpha >  tol_alpha;
% idx(:,:,6)  = abs(S - 1) >  tol_S & abs(T) <= tol_T & abs(u) >  tol_u &  alpha <= tol_alpha;
% idx(:,:,7)  = abs(S - 1) >  tol_S & abs(T) <= tol_T & abs(u) <= tol_u &  alpha >  tol_alpha;
% idx(:,:,8)  = abs(S - 1) >  tol_S & abs(T) <= tol_T & abs(u) <= tol_u &  alpha <= tol_alpha;
% idx(:,:,9)  = abs(S - 1) <= tol_S & abs(T) >  tol_T & abs(u) >  tol_u &  alpha >  tol_alpha;
% idx(:,:,10) = abs(S - 1) <= tol_S & abs(T) >  tol_T & abs(u) >  tol_u &  alpha <= tol_alpha;
% idx(:,:,11) = abs(S - 1) <= tol_S & abs(T) >  tol_T & abs(u) <= tol_u &  alpha >  tol_alpha;
% idx(:,:,12) = abs(S - 1) <= tol_S & abs(T) >  tol_T & abs(u) <= tol_u &  alpha <= tol_alpha;
% idx(:,:,13) = abs(S - 1) <= tol_S & abs(T) <= tol_T & abs(u) >  tol_u &  alpha >  tol_alpha;
% idx(:,:,14) = abs(S - 1) <= tol_S & abs(T) <= tol_T & abs(u) >  tol_u &  alpha <= tol_alpha;
% idx(:,:,15) = abs(S - 1) <= tol_S & abs(T) <= tol_T & abs(u) <= tol_u &  alpha >  tol_alpha;
% idx(:,:,16) = abs(S - 1) <= tol_S & abs(T) <= tol_T & abs(u) <= tol_u &  alpha <= tol_alpha;
% idx(:,:,17) = abs(S - 1) >  tol_S & abs(T) >  tol_T & abs(u) >  tol_u & abs((u - alpha)./alpha) <=  tol_ualpha;
% idx(:,:,18) = abs(S - 1) >  tol_S & abs(T) <= tol_T & abs(u) >  tol_u & abs((u - alpha)./alpha) <=  tol_ualpha;
% idx(:,:,19) = abs(S - 1) <= tol_S & abs(T) >  tol_T & abs(u) >  tol_u & abs((u - alpha)./alpha) <=  tol_ualpha;
% idx(:,:,20) = abs(S - 1) <= tol_S & abs(T) <= tol_T & abs(u) >  tol_u & abs((u - alpha)./alpha) <=  tol_ualpha;

% idx(:,:,1) = abs(u) > tol_u & alpha > tol_alpha;
% idx(:,:,2) = abs(u) > tol_u & alpha <= tol_alpha;
% idx(:,:,3)  = abs(u) <= tol_u &  alpha >  tol_alpha;
% idx(:,:,4)  = abs(u) <= tol_u &  alpha <= tol_alpha;
% idx(:,:,17) = abs(u) > tol_u & abs((u - alpha)./alpha) <= tol_ualpha;

idx(:,:,1) = alpha > tol_alpha;
idx(:,:,2) = alpha <= tol_alpha;
% idx(:,:,17) = abs((u - alpha)./alpha) <= tol_ualpha;

K0 = fcnK0(S, T, u, alpha, F_lim, tol, idx);
K1 = fcnK1(S, T, u, alpha, F_lim, tol, idx);
K2 = fcnK2(S, T, u, alpha, F_lim, tol, idx);
K3 = fcnK3(S, T, u, alpha, F_lim, tol, idx);
K4 = fcnK4(S, T, u, alpha, F_lim, tol, idx);
K5 = fcnK5(S, T, u, alpha, F_lim, tol, idx);
K6 = fcnK6(S, T, u, alpha, F_lim, tol, idx);
K7 = fcnK7(S, T, u, alpha, F_lim, tol, idx);

tol = 1e-2;
len = size(S,1);
for i = 1:len
    Ka = @(F) Kaf(F, S(i,1), T(i,1), u(i,1), alpha(i,1));
    if sign(F_lim(i,1,1)) ~= sign(F_lim(i,1,2))
        K0n(i,1) = integral(Ka, F_lim(i,1,1), -tol) + integral(Ka, tol, F_lim(i,1,2));
    else
        K0n(i,1) = integral(Ka, F_lim(i,1,1), F_lim(i,1,2));
    end
end

figure(1)
clf(1) 
scatter3(C(:,1), N(:,1), K0(:,1), 'ko','linewidth',2);
hold on
scatter3(C(:,1), N(:,1), K0n(:,1), 'sr','linewidth',2);
hold off
grid minor
box on
axis tight
legend('Analytical','Numerica','Location','NorthEast')
xlabel('C');
ylabel('N');
zlabel('Int Value');

function Ka_out = Kaf(F, S, T, u, alpha)

alpha = repmat(alpha, 1, size(F,2));
u = repmat(u, 1, size(F,2));

denom = (((S.*(F.^2) + (F.*T) + u).^(3/2)).*(((F.^2) + alpha).^2));

Ka_out = 1./denom;

end


