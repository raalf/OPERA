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

W = sqrt(alpha.^2 + (b.^2 - 2.*a).*alpha + a.^2);
K = W + a - alpha;
L = K./(b.^2);
G = a./(b.^2);
D = alpha./(b.^2);

p = 1 - ((2.*a)./K) + ((2.*alpha)./K);
q = ((6.*D - 6.*G + 2.*L + 1).*K + (4.*G - 1).*a + (4.*D - 8.*G + 1).*alpha)./...
    ((2.*D - 2.*G + 2.*L - 1).*K + a - alpha);

mu = -M.*K + 2.*M.*a - 2.*M.*alpha - N.*b;
lambda = K.*M - N.*b;

t1 = (-F1.*b + sqrt(alpha.*b.^2 + a.^2 - 2.*a.*alpha + alpha.^2) - a + alpha)./...
     (F1.*b + sqrt(alpha.*b.^2 + a.^2 - 2.*a.*alpha + alpha.^2) + a - alpha);

t2 = (-F2.*b + sqrt(alpha.*b.^2 + a.^2 - 2.*a.*alpha + alpha.^2) - a + alpha)./...
     (F2.*b + sqrt(alpha.*b.^2 + a.^2 - 2.*a.*alpha + alpha.^2) + a - alpha); 
 
t1(abs(t1) > 1e+10,1) = sign(t1(abs(t1) > 1e+10)).*1e+10;
t2(abs(t2) > 1e+10,1) = sign(t2(abs(t2) > 1e+10)).*1e+10;

I = lambda.*-(1./K).*fcnH_3(p, q, t1, t2) + mu.*(1./K).*fcnH_4(p, q, t1, t2);

end