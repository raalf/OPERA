function I = fcnK5G(S, T, u, alpha, F, tol)
t1 = F .^ 2;t3 = sqrt(alpha);t4 = 0.1e1 ./ t3;t6 = atan(t4 .* F);t12 = t1 .* S;t13 = sqrt(t12);t22 = -t4 ./ (t1 + alpha) ./ t13 ./ t12 .* (t3 .* F - alpha .* t6 - t1 .* t6) .* t1 .* F ./ 0.2e1;
I = t22(:,:,2) - t22(:,:,1);

end
