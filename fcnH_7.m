function I = fcnH_7(M, N, alpha, F1, F2)
% Solves an integral of the form,
%
%  F2
%   /
%  |                N + M F
%  |  arcsinh(------------------) dF
% /                  2
%  F1          sqrt(F  + alpha) 
%

%%
len = size(N,1);
I = nan(len,1);
tol = 1e-5;
for q = 1:len
    num_eqn = @(F) asinh((N(q) + M(q).*F)./sqrt((F.^2 + alpha(q))));
    if sign(F1(q)) ~= sign(F2(q)) && abs(alpha(q)) < 1e-5
        I(q) = integral(num_eqn, F1(q), sign(F1(q)).*tol) + integral(num_eqn, sign(F2(q)).*tol, F2(q));
    else
        I(q) = integral(num_eqn, F1(q), F2(q));
    end
end

end