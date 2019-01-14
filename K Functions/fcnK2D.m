function I = fcnK2D(S, T, u, alpha, F, tol)
t1 = (F .* S);t3 = (S .^ 2);t5 = (F .^ 2);t12 = (T .^ 2);t19 = t12 .^ 2;t24 = F .* T + S .* t5;t25 = sqrt(t24);t31 = -0.2e1 ./ 0.5e1 ./ t25 ./ t24 ./ t19 ./ F .* (16 .* t5 .* F .* t3 .* S + 8 .* T .* t5 .* t3 + t12 .* T - 2 .* t12 .* t1) .* (t1 + T);
I = t31(:,:,2) - t31(:,:,1);

end
