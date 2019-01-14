function I = fcnK2L(S, T, u, alpha, F, tol)
t2 = (F .^ 2);t7 = (T .^ 2);t14 = t7 .^ 2;t18 = F .* T + t2;t19 = sqrt(t18);t25 = -0.2e1 ./ 0.5e1 ./ t19 ./ t18 ./ t14 ./ F .* (16 .* t2 .* F - 2 .* t7 .* F + 8 .* T .* t2 + t7 .* T) .* (F + T);
I = t25(:,:,2) - t25(:,:,1);

end
