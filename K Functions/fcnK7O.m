function I = fcnK7O(S, T, u, alpha, F, tol)
t1 = F .^ 2;t2 = t1 .* F;t3 = sqrt(alpha);t6 = 0.1e1 ./ t3;t8 = atan(t6 .* F);t15 = alpha .^ 2;t20 = sqrt(t1);t28 = t6 ./ (t1 + alpha) ./ t20 ./ t1 .* (0.3e1 .* t3 .* alpha .* F - 0.3e1 .* alpha .* t1 .* t8 - 0.3e1 .* t15 .* t8 + 0.2e1 .* t2 .* t3) .* t2 ./ 0.2e1;
I = t28(:,:,2) - t28(:,:,1);

end
