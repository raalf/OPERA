function I = fcnK4D(S, T, u, alpha, F, tol)
t1 = (F .* S);t6 = T .^ 2;t9 = F .^ 2;t12 = F .* T + S .* t9;t13 = sqrt(t12);t19 = -0.2e1 ./ t13 ./ t12 ./ t6 .* (2 .* t1 + T) .* (t1 + T) .* F;
I = t19(:,:,2) - t19(:,:,1);

end
