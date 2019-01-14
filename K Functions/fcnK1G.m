function I = fcnK1G(S, T, u, alpha, F, tol)
t1 = F .^ 2;t2 = sqrt(alpha);t5 = atan(0.1e1 ./ t2 .* F);t18 = t1 .* S;t19 = sqrt(t18);t22 = alpha .^ 2;t31 = -0.1e1 ./ (t1 + alpha) ./ t2 ./ t22 ./ t19 ./ t18 .* (0.3e1 .* alpha .* F .* t5 + 0.3e1 .* t1 .* F .* t5 + 0.2e1 .* alpha .* t2 + 0.3e1 .* t2 .* t1) .* t1 ./ 0.2e1;
I = t31(:,:,2) - t31(:,:,1);

end
