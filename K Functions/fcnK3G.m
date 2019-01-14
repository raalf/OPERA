function I = fcnK3G(S, T, u, alpha, F, tol)
t1 = F .^ 2;t3 = sqrt(alpha);t8 = atan(0.1e1 ./ t3 .* F);t11 = alpha .^ 2;t15 = t1 .* S;t16 = sqrt(t15);t26 = 0.1e1 ./ (t1 + alpha) ./ t3 ./ t11 ./ t16 ./ t15 .* (t3 .* alpha .* F + alpha .* t1 .* t8 + t11 .* t8) .* t1 .* F ./ 0.2e1;
I = t26(:,:,2) - t26(:,:,1);

end
