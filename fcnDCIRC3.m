function [circ_int] = fcnDCIRC3(xi_len, eta_len, valWNELE)
% For finding integrated circulation

c3eta = eta_len;
a2 = 0.5.*eta_len.^2;
a1 = (1.6).*eta_len.^3;

c3xi = xi_len;
b2 = 0.5.*xi_len.^2;
b1 = (1.6).*xi_len.^3;

gamma1 = [a1, a2, zeros(size(eta_len)), zeros(size(eta_len)), c3eta];
gamma2 = [zeros(size(eta_len)), zeros(size(eta_len)), b1, b2, c3xi];
circ_int = [fcnCREATEDSECT(sparse(length(eta_len), valWNELE*5), length(eta_len), 5, [1:valWNELE]', [], gamma1, []); fcnCREATEDSECT(sparse(length(eta_len), valWNELE*5), length(eta_len), 5, [1:valWNELE]', [], gamma2, [])];



end

