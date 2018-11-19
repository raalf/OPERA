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

len = size(N,1);

b = T./S;
a = u./S;

tau = (2.*sqrt((a - alpha).^2 + alpha.*(b.^2)))./(-b);
lambda = -M.*((a - alpha)./b) + ((M.*tau)./2) + N;
mu = -M.*((a - alpha)./b) - ((M.*tau)./2) + N;
phi = (((a - alpha)./b) + (-tau./2)).^2 - alpha;
psi = phi - sqrt((a - alpha).^2 + (alpha.*(b.^2)));
p = ((((a - alpha)./b) + (tau./2)).^2 + alpha)./phi;
q = ((((a - alpha)./b) + (tau./2)).^2 + sqrt((a - alpha).^2 + (alpha.*(b.^2))) + alpha)./psi;

% idx = (abs(t) > 1e10) | (N == 0) | (q < -(t.^2) | q.^2 < -(t.^2)) | ((q > p) & (p < 0));  % Correction
% tau(idx) = nan; % Correction
% lambda(idx) = nan; % Correction
% mu(idx) = nan;  % Correction
% phi(idx) = nan; % Correction
% psi(idx) = nan; % Correction
% p(idx) = nan; % Correction
% q(idx) = nan; % Correction
% t(idx) = nan; % Correction

tdt = nan(len,2);
dt = nan(len,2);

%% At F1
t(:,1) = (((a - alpha)./b) + (tau./2) + F1)./(((a - alpha)./-b) + (tau./2) - F1);
t(:,2) = (((a - alpha)./b) + (tau./2) + F2)./(((a - alpha)./-b) + (tau./2) - F2);

idx = p > q;
tdt(idx,:) = (1./(sqrt(p(idx) - q(idx)))).*atan2(real(sqrt(t(idx,:).^2 + q(idx))), sqrt(p(idx) - q(idx)));
dt(idx,:) = (1./(2.*sqrt(p(idx)).*sqrt(q(idx) - p(idx)))).*log((sqrt(p(idx)).*(sqrt((t(idx,:).^2) + q(idx))) + t(idx,:).*(sqrt(p(idx) - q(idx))))./(sqrt(p(idx)).*(sqrt((t(idx,:).^2) + q(idx))) - t(idx,:).*(sqrt(p(idx) - q(idx)))));

idx = q > p;
tdt(idx) = (-1./(2.*sqrt(q(idx) - p(idx)))).*log((sqrt(q(idx) - p(idx)) + sqrt(t(idx).^2 + q(idx)))./(sqrt(q(idx) - p(idx)) - sqrt(t(idx).^2 + q(idx))));
dt(idx) = (1./(sqrt(p(idx)).*sqrt(q(idx) - p(idx)))).*atan2(real(t(idx).*sqrt(q(idx) - p(idx))), real(sqrt(p(idx)).*sqrt(t(idx).^2 + q(idx).^2)));

I1 = sign(t(:,1) + 1).*(tau./(phi.*sqrt(S.*psi))).*(lambda.*tdt(:,1) + mu.*dt(:,1));
I2 = sign(t(:,2) + 1).*(tau./(phi.*sqrt(S.*psi))).*(lambda.*tdt(:,2) + mu.*dt(:,2));

%%
I = I2 - I1;
I = real(I);
end