function I = fcnK6B(S, T, u, alpha, F, tol)
t1 = sqrt(S);t2 = t1 .* S;t4 = F .^ 2;t6 = F .* T;t8 = sqrt(t4 .* S + t6 + u);t17 = S .* u;t18 = T .^ 2;t22 = log(0.2e1);t30 = log(0.1e1 ./ t1 .* (0.2e1 .* S .* F + 0.2e1 .* t1 .* t8 + T));t39 = 0.4e1 ./ (0.4e1 .* t17 - t18) .* (-u .* t2 .* F + (t6 + u) .* t1 .* T ./ 0.2e1 + (-t22 + t30) .* (t17 - t18 ./ 0.4e1) .* t8) ./ t8 ./ t2;
I = t39(:,:,2) - t39(:,:,1);

end
