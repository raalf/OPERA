function I = fcnK2O(S, T, u, alpha, F, tol)
t1 = F .^ 2;t3 = log(F);t6 = t1 + alpha;t7 = log(t6);t14 = sqrt(t1);t19 = alpha .^ 2;t23 = 0.1e1 ./ t19 ./ t6 ./ t14 .* (0.2e1 .* alpha .* t3 - alpha .* t7 + 0.2e1 .* t1 .* t3 - t1 .* t7 + alpha) .* F ./ 0.2e1;
I = t23(:,:,2) - t23(:,:,1);

end
