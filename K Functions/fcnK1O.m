function I = fcnK1O(S, T, u, alpha, F, tol)
t1 = F .^ 2;t2 = sqrt(alpha);t5 = atan(0.1e1 ./ t2 .* F);t18 = sqrt(t1);t21 = alpha .^ 2;t30 = -0.1e1 ./ (t1 + alpha) ./ t21 ./ t2 ./ t18 .* (0.3e1 .* alpha .* F .* t5 + 0.3e1 .* t1 .* F .* t5 + 0.2e1 .* alpha .* t2 + 0.3e1 .* t2 .* t1) ./ 0.2e1;
I = t30(:,:,2) - t30(:,:,1);

end
