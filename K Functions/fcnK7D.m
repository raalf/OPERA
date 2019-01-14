function I = fcnK7D(S, T, u, alpha, F, tol)
t1 = F .^ 2;t7 = sqrt(F .* T + S .* t1);t8 = 0.1e1 ./ t7;t10 = S .^ 2;t16 = sqrt(S);t26 = log(0.1e1 ./ t16 .* (T ./ 0.2e1 + F .* S) + t7);t29 = t8 ./ S .* t1 + 0.3e1 .* F .* t8 ./ t10 .* T - 0.3e1 ./ 0.2e1 .* t26 ./ t16 ./ t10 .* T;
I = t29(:,:,2) - t29(:,:,1);

end
