function I = fcnK1L(S, T, u, alpha, F, tol)
t2 = (F .^ 2);t3 = (t2 .^ 2);t8 = (T .^ 2);t14 = (t8 .^ 2);t23 = F .* T + t2;t24 = sqrt(t23);t29 = 0.2e1 ./ 0.35e2 ./ t24 ./ t23 ./ t14 ./ T ./ t2 .* (64 .* T .* t2 .* F + 8 .* t8 .* T .* F - 16 .* t8 .* t2 - 5 .* t14 + 128 .* t3) .* (F + T);
I = t29(:,:,2) - t29(:,:,1);

end
