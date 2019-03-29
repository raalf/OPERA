clc
clear
warning off
tol = 1e-10;

len = 100;
%%
F_lim(:,:,1) = zeros(len, 1) + 2;
F_lim(:,:,2) = zeros(len, 1) + 1;

% F_lim(:,:,1) = zeros(len, 1) - 2;
% F_lim(:,:,2) = zeros(len, 1) - 0.000001;
 
S = linspace(2,5,len)';
T = zeros(len, 1) + 1;
u = zeros(len, 1) + 1;
alpha = zeros(len, 1) ;

%%
F_lim = repmat(F_lim, 1, 1, 1);
S = repmat(S, 1, 2);
T = repmat(T, 1, 2);
u = repmat(u, 1, 2);

tol_F = 1e-10;
tol_S = 1e-4;
tol_T = 1e-4;
tol_u = 1e-10;
tol_alpha = 1e-10;
tol_ualpha = 1e-3;

% Correcting zero F limits if we are in the plane of the element
idx = abs(F_lim(:,:,1)) < tol_F & alpha(:,1) <= 1e-3;
sgn = sign(F_lim(idx,:,1));
sgn(sign(F_lim(idx,:,1)) == 0) = sign(F_lim(sign(F_lim(idx,:,1)) == 0,:,2));
% F_lim(idx,:,1) = sgn.*tol_F;
% alpha(idx) = 1e-1;

idx = abs(F_lim(:,:,2)) < tol_F & alpha <= 1e-3;
sgn = sign(F_lim(idx,:,2));
sgn(sign(F_lim(idx,:,2)) == 0) = sign(F_lim(sign(F_lim(idx,:,2)) == 0,:,1));
% F_lim(idx,:,2) = sgn.*tol_F;
% alpha(idx,:) = 1e-1;

F_lim = repmat(F_lim,1,2,1);
alpha = repmat(alpha,1,2,1);

% Identifying which cases are special
idx = [];
idx(:,:,1)  = abs(S - 1) >  tol_S & abs(T) >  tol_T & abs(u) >  tol_u &  alpha >  tol_alpha;
idx(:,:,2)  = abs(S - 1) >  tol_S & abs(T) >  tol_T & abs(u) >  tol_u &  alpha <= tol_alpha;
idx(:,:,3)  = abs(S - 1) >  tol_S & abs(T) >  tol_T & abs(u) <= tol_u &  alpha >  tol_alpha;
idx(:,:,4)  = abs(S - 1) >  tol_S & abs(T) >  tol_T & abs(u) <= tol_u &  alpha <= tol_alpha;
idx(:,:,5)  = abs(S - 1) >  tol_S & abs(T) <= tol_T & abs(u) >  tol_u &  alpha >  tol_alpha;
idx(:,:,6)  = abs(S - 1) >  tol_S & abs(T) <= tol_T & abs(u) >  tol_u &  alpha <= tol_alpha;
idx(:,:,7)  = abs(S - 1) >  tol_S & abs(T) <= tol_T & abs(u) <= tol_u &  alpha >  tol_alpha;
idx(:,:,8)  = abs(S - 1) >  tol_S & abs(T) <= tol_T & abs(u) <= tol_u &  alpha <= tol_alpha;
idx(:,:,9)  = abs(S - 1) <= tol_S & abs(T) >  tol_T & abs(u) >  tol_u &  alpha >  tol_alpha;
idx(:,:,10) = abs(S - 1) <= tol_S & abs(T) >  tol_T & abs(u) >  tol_u &  alpha <= tol_alpha;
idx(:,:,11) = abs(S - 1) <= tol_S & abs(T) >  tol_T & abs(u) <= tol_u &  alpha >  tol_alpha;
idx(:,:,12) = abs(S - 1) <= tol_S & abs(T) >  tol_T & abs(u) <= tol_u &  alpha <= tol_alpha;
idx(:,:,13) = abs(S - 1) <= tol_S & abs(T) <= tol_T & abs(u) >  tol_u &  alpha >  tol_alpha;
idx(:,:,14) = abs(S - 1) <= tol_S & abs(T) <= tol_T & abs(u) >  tol_u &  alpha <= tol_alpha;
idx(:,:,15) = abs(S - 1) <= tol_S & abs(T) <= tol_T & abs(u) <= tol_u &  alpha >  tol_alpha;
idx(:,:,16) = abs(S - 1) <= tol_S & abs(T) <= tol_T & abs(u) <= tol_u &  alpha <= tol_alpha;
idx(:,:,17) = abs(S - 1) >  tol_S & abs(T) >  tol_T & abs(u) >  tol_u & abs(u - alpha) <=  tol_u;
idx(:,:,18) = abs(S - 1) >  tol_S & abs(T) <= tol_T & abs(u) >  tol_u & abs(u - alpha) <=  tol_u;
idx(:,:,19) = abs(S - 1) <= tol_S & abs(T) >  tol_T & abs(u) >  tol_u & abs(u - alpha) <=  tol_u;
idx(:,:,20) = abs(S - 1) <= tol_S & abs(T) <= tol_T & abs(u) >  tol_u & abs(u - alpha) <=  tol_u;
% idx(:,:,17) = abs(S - 1) >  tol_S & abs(T) >  tol_T & abs(u) >  tol_u & abs((u - alpha)./alpha) <=  tol_ualpha;
% idx(:,:,18) = abs(S - 1) >  tol_S & abs(T) <= tol_T & abs(u) >  tol_u & abs((u - alpha)./alpha) <=  tol_ualpha;
% idx(:,:,19) = abs(S - 1) <= tol_S & abs(T) >  tol_T & abs(u) >  tol_u & abs((u - alpha)./alpha) <=  tol_ualpha;
% idx(:,:,20) = abs(S - 1) <= tol_S & abs(T) <= tol_T & abs(u) >  tol_u & abs((u - alpha)./alpha) <=  tol_ualpha;

K0 = fcnK0(S, T, u, alpha, F_lim, tol, idx);
% K1 = fcnK1(S, T, u, alpha, F_lim, tol, idx);
% K2 = fcnK2(S, T, u, alpha, F_lim, tol, idx);
% K3 = fcnK3(S, T, u, alpha, F_lim, tol, idx);
% K4 = fcnK4(S, T, u, alpha, F_lim, tol, idx);
% K5 = fcnK5(S, T, u, alpha, F_lim, tol, idx);
% K6 = fcnK6(S, T, u, alpha, F_lim, tol, idx);
% K7 = fcnK7(S, T, u, alpha, F_lim, tol, idx);


tol = 1e-10;
for i = 1:len
%     i
%     syms F
%     Ka = 1./(((S(i,1).*(F.^2) + (F.*T(i,1)) + u(i,1)).^(3/2)).*(((F.^2) + alpha(i,1)).^2));
%     K0n(i,1) = vpaintegral(Ka, F, F_lim(i,1,1), F_lim(i,1,2));

    Ka = @(F) Kaf(F, S(i,1), T(i,1), u(i,1), alpha(i,1));
    if sign(F_lim(i,1,1)) ~= sign(F_lim(i,1,2))
        K0n(i,1) = integral(Ka, F_lim(i,1,1), 0) + integral(Ka, 0, F_lim(i,1,2));
    else
        K0n(i,1) = integral(Ka, F_lim(i,1,1), F_lim(i,1,2));
    end
end

figure(1)
clf(1)
plot(1:len, K0(:,1), '-k');
hold on
plot(1:len, K0n, '--b');
hold off
grid minor
box on
axis tight

function Ka_out = Kaf(F, S, T, u, alpha)

alpha = repmat(alpha, 1, size(F,2));
u = repmat(u, 1, size(F,2));

denom = (((S.*(F.^2) + (F.*T) + u).^(3/2)).*(((F.^2) + alpha).^2));

alpha(abs(denom) < 1e-5) = alpha(abs(denom) < 1e-5) + 1e-1;
u(abs(denom) < 1e-5) = u(abs(denom) < 1e-5) + 2e-1;
denom = (((S.*(F.^2) + (F.*T) + u).^(3/2)).*(((F.^2) + alpha).^2));

Ka_out = 1./denom;

end


