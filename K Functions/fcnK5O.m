function I = fcnK5O(S, T, u, alpha, F, tol)
t1 = F .^ 2;t3 = sqrt(alpha);t4 = 0.1e1 ./ t3;t6 = atan(t4 .* F);t12 = sqrt(t1);t21 = -t4 ./ (t1 + alpha) ./ t12 .* (t3 .* F - alpha .* t6 - t6 .* t1) .* F ./ 0.2e1;
I = t21(:,:,2) - t21(:,:,1);

end
