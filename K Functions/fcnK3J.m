function I = fcnK3J(S, T, u, alpha, F, tol)
t1 = 0.1e1 ./ u;t2 = F .^ 2;t3 = F .* T;t5 = sqrt(t2 + t3 + u);t6 = 0.1e1 ./ t5;t11 = T .^ 2;t18 = sqrt(u);t27 = log(0.1e1 ./ F .* (0.2e1 .* t5 .* t18 + t3 + 0.2e1 .* u));t29 = t6 .* t1 - t6 ./ (-t11 + 0.4e1 .* u) .* (0.2e1 .* F + T) .* t1 .* T - t27 ./ t18 ./ u;
I = t29(:,:,2) - t29(:,:,1);

end
