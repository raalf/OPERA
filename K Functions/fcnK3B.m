function I = fcnK3B(S, T, u, alpha, F, tol)
t1 = sqrt(u);t2 = t1 .* u;t10 = F .* T;t11 = F .^ 2;t14 = sqrt(t11 .* S + t10 + u);t20 = log(0.1e1 ./ F .* (0.2e1 .* t14 .* t1 + t10 + 0.2e1 .* u));t24 = T .^ 2;t25 = 0.4e1 .* S .* u - t24;t36 = -0.4e1 ./ t25 ./ t14 ./ t2 .* (-t2 .* S + (S .* F + T) .* t1 .* T ./ 0.2e1 + t25 .* t14 .* t20 ./ 0.4e1);
I = t36(:,:,2) - t36(:,:,1);

end
