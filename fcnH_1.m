function I = fcnH_1(M, N, S, T, u, alpha, F1, F2)
% Solves an integral of the form,
%
%  F2
%   /
%  |               M F + N
%  |  --------------------------------- dF
% /     2                  2
%  F1 (F  + alpha) sqrt(S F  + T F + u)
%
% Using the method of David Yeh
% https://books.google.ca/books?id=F7jiBQAAQBAJ&printsec=frontcover#v=onepage&q&f=false

%%
% len = size(N,1);
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
% dt(idx,:) = (1./(2.*sqrt(p(idx)).*sqrt(q(idx) - p(idx)))).*log(abs((sqrt(p(idx)).*(sqrt((t(idx,:).^2) + q(idx))) + t(idx,:).*(sqrt(p(idx) - q(idx))))./(sqrt(p(idx)).*(sqrt((t(idx,:).^2) + q(idx))) - t(idx,:).*(sqrt(p(idx) - q(idx))))));
% 
% idx = q > p;
% tdt(idx,:) = (-1./(2.*sqrt(q(idx) - p(idx)))).*log(abs((sqrt(q(idx) - p(idx)) + sqrt(t(idx,:).^2 + q(idx)))./(sqrt(q(idx) - p(idx)) - sqrt(t(idx,:).^2 + q(idx)))));
% dt(idx,:) = (1./(sqrt(p(idx)).*sqrt(q(idx) - p(idx)))).*atan2(real(t(idx,:).*sqrt(q(idx) - p(idx))), real(sqrt(p(idx)).*sqrt(t(idx,:).^2 + q(idx).^2)));
% 
% I1 = sign(t(:,1) + 1).*(tau./(phi.*sqrt(S.*psi))).*(lambda.*tdt(:,1) + mu.*dt(:,1));
% I2 = sign(t(:,2) + 1).*(tau./(phi.*sqrt(S.*psi))).*(lambda.*tdt(:,2) + mu.*dt(:,2));
% I = I2 - I1;

%%
len = size(N,1);
I = nan(len,1);
tol = 8e-2;
for q = 1:len
    num_eqn = @(F) (N(q) + M(q).*F)./((F.^2 + alpha(q)).*sqrt(S(q).*F.^2 + T(q).*F + u(q)));
    if sign(F1(q)) ~= sign(F2(q)) && abs(alpha(q)) < 1e-5
        I(q) = integral(num_eqn, F1(q), sign(F1(q)).*tol) + integral(num_eqn, sign(F2(q)).*tol, F2(q));
    else 
        I(q) = integral(num_eqn, F1(q), F2(q));
    end
end

end