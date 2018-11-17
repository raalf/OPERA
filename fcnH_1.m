function I = fcnH_1(C, N, alpha, F)
% Solves an integral of the form, 
%
%  F2
%   /
%  |               N + C F
%  |  --------------------------------- dF
% /     2                  2
%  F1 (F  + alpha) sqrt(S F  + T F + u)
%
% Using the method of David Yeh


S = C.^2 + 1;
T = 2.*C.*N;
u = N.^2 + alpha;

b = T./S;
a = u./S;
tau = (2.*sqrt((a - alpha).^2 + alpha.*(b.^2)))./(-b);
lambda = -C.*((a - alpha)./b) + ((C.*tau)./2) + N;
mu = -C.*((a - alpha)./b) - ((C.*tau)./2) + N;
phi = (((a - alpha)./b) + (-tau./2)).^2 - alpha;
psi = phi - sqrt((a - alpha).^2 + (alpha.*(b.^2)));
p = ((((a - alpha)./b) + (tau./2)).^2 + alpha)./phi;
q = ((((a - alpha)./b) + (tau./2)).^2 + sqrt((a - alpha).^2 + (alpha.*(b.^2))) + alpha)./psi;
t = (((a - alpha)./b) + (tau./2) + F)./(((a - alpha)./-b) + (tau./2) - F);

idx = (abs(t) > 1e10) | (N == 0) | (q < -(t.^2) | q.^2 < -(t.^2)) | ((q > p) & (p < 0));  % Correction
tau(idx) = nan; % Correction
lambda(idx) = nan; % Correction
mu(idx) = nan;  % Correction
phi(idx) = nan; % Correction
psi(idx) = nan; % Correction
p(idx) = nan; % Correction
q(idx) = nan; % Correction
t(idx) = nan; % Correction

len = size(N,1);
tdt = nan(len,1);
dt = nan(len,1);

idx = p > q;
tdt(idx) = (1./(sqrt(p(idx) - q(idx)))).*atan2(sqrt(t(idx).^2 + q(idx)), sqrt(p(idx) - q(idx)));
dt(idx) = (1./(2.*sqrt(p(idx)).*sqrt(q(idx) - p(idx)))).*log((sqrt(p(idx)).*(sqrt((t(idx).^2) + q(idx))) + t(idx).*(sqrt(p(idx) - q(idx))))./(sqrt(p(idx)).*(sqrt((t(idx).^2) + q(idx))) - t(idx).*(sqrt(p(idx) - q(idx)))));

idx = q > p;
tdt(idx) = (-1./(2.*sqrt(q(idx) - p(idx)))).*log((sqrt(q(idx) - p(idx)) + sqrt(t(idx).^2 + q(idx)))./(sqrt(q(idx) - p(idx)) - sqrt(t(idx).^2 + q(idx))));
dt(idx) = (1./(sqrt(p(idx)).*sqrt(q(idx) - p(idx)))).*atan2(t(idx).*sqrt(q(idx) - p(idx)), sqrt(p(idx)).*sqrt(t(idx).^2 + q(idx).^2));

I = sign(t + 1).*(tau./(phi.*sqrt(S.*psi))).*(lambda.*tdt + mu.*dt);

I = real(I); % Correction

end