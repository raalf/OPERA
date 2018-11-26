function I = fcnH_1(M, N, S, T, u, alpha, F1, F2)
% Solves an integral of the form,
%
%  F2
%   /
%  |               N + M F
%  |  --------------------------------- dF
% /     2                  2
%  F1 (F  + alpha) sqrt(S F  + T F + u)
%
% Using the method of David Yeh
% https://books.google.ca/books?id=F7jiBQAAQBAJ&printsec=frontcover#v=onepage&q&f=false

%%
% len = size(N,1);
% 
% 
% b = T./S;
% a = u./S;
% tau = (2.*sqrt((a - alpha).^2 + alpha.*b.^2))./-b;
% lambda = -M.*((a - alpha)./b) + ((M.*tau)./2) + N;
% mu = -M.*((a - alpha)./b) - ((M.*tau)./2) + N;
% phi = (((a - alpha)./b) + (-tau./2)).^2 - alpha;
% psi = phi - sqrt((a - alpha).^2 + alpha.*(b.^2));
% p = ((((a - alpha)./b) + (tau./2)).^2 + alpha)./phi;
% q = ((((a - alpha)./b) + (tau./2)).^2 + sqrt((a - alpha).^2 + alpha.*(b.^2)) + alpha)./psi;
% t = (((a - alpha)./b) + (tau./2) + [F1 F2])./(((a - alpha)./-b) + (tau./2) - [F1 F2]);
% 
% tdt = nan(len,2);
% dt = nan(len,2);
% 
% idx = p > q;
% tdt(idx,:) = (1./(sqrt(p(idx) - q(idx)))).*atan2(real(sqrt(t(idx,:).^2 + q(idx))), sqrt(p(idx) - q(idx)));
% dt(idx,:) = (1./(2.*sqrt(p(idx)).*sqrt(q(idx) - p(idx)))).*log((sqrt(p(idx)).*(sqrt((t(idx,:).^2) + q(idx))) + t(idx,:).*(sqrt(p(idx) - q(idx))))./(sqrt(p(idx)).*(sqrt((t(idx,:).^2) + q(idx))) - t(idx,:).*(sqrt(p(idx) - q(idx)))));
% 
% idx = q > p;
% tdt(idx,:) = (-1./(2.*sqrt(q(idx) - p(idx)))).*log((sqrt(q(idx) - p(idx)) + sqrt(t(idx,:).^2 + q(idx)))./(sqrt(q(idx) - p(idx)) - sqrt(t(idx,:).^2 + q(idx))));
% dt(idx,:) = (1./(sqrt(p(idx)).*sqrt(q(idx) - p(idx)))).*atan2(real(t(idx,:).*sqrt(q(idx) - p(idx))), real(sqrt(p(idx)).*sqrt(t(idx,:).^2 + q(idx).^2)));
% 
% I1 = sign(t(:,1) + 1).*(tau./(phi.*sqrt(S.*psi))).*(lambda.*tdt(:,1) + mu.*dt(:,1));
% I2 = sign(t(:,2) + 1).*(tau./(phi.*sqrt(S.*psi))).*(lambda.*tdt(:,2) + mu.*dt(:,2));
% I = I2 - I1;

%%
len = size(N,1);
I = nan(len,1);
tol = 8e-4;
for q = 1:len
    num_eqn = @(F) (N(q) + M(q).*F)./((F.^2 + alpha(q)).*sqrt(S(q).*F.^2 + T(q).*F + u(q)));
    if sign(F1(q)) ~= sign(F2(q)) && abs(alpha(q)) < 1e-5
        I(q) = integral(num_eqn, F1(q), sign(F1(q)).*tol) + integral(num_eqn, sign(F2(q)).*tol, F2(q));
    else 
        I(q) = integral(num_eqn, F1(q), F2(q));
    end
end

% W = sqrt(alpha.^2 + (b.^2 - 2.*a).*alpha + a.^2);
% K = W + a - alpha;
% L = K./(b.^2);
% G = a./(b.^2);
% D = alpha./(b.^2);
% r = S.*((2.*D - 2.*G + 2.*L - 1).*K + a - alpha);
%
% p = 1 - ((2.*a)./K) + ((2.*alpha)./K);
% q = ((6.*D - 6.*G + 2.*L + 1).*K + (4.*G - 1).*a + (4.*D - 8.*G + 1).*alpha)./...
%     ((2.*D - 2.*G + 2.*L - 1).*K + a - alpha);
%
% mu = (-M.*K + 2.*M.*a - 2.*M.*alpha - N.*b)./(K.*sqrt(r));
% lambda = (K.*M - N.*b)./(K.*sqrt(r));
%
% t1 = (-F1.*b + sqrt(alpha.*b.^2 + a.^2 - 2.*a.*alpha + alpha.^2) - a + alpha)./...
%      (F1.*b + sqrt(alpha.*b.^2 + a.^2 - 2.*a.*alpha + alpha.^2) + a - alpha);
%
% t2 = (-F2.*b + sqrt(alpha.*b.^2 + a.^2 - 2.*a.*alpha + alpha.^2) - a + alpha)./...
%      (F2.*b + sqrt(alpha.*b.^2 + a.^2 - 2.*a.*alpha + alpha.^2) + a - alpha);
%
% t1(abs(t1) > 1e+10,1) = sign(t1(abs(t1) > 1e+10)).*1e+10;
% t2(abs(t2) > 1e+10,1) = sign(t2(abs(t2) > 1e+10)).*1e+10;
%
% I = -lambda.*fcnH_3(p, q, t1, t2) + mu.*fcnH_4(p, q, t1, t2);

end