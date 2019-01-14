function I = fcnK6G(S, T, u, alpha, F, tol)
t1 = F .^ 2;t3 = t1 + alpha;t4 = log(t3);t10 = t1 .* S;t11 = sqrt(t10);t16 = 0.1e1 ./ t11 ./ t10 ./ t3 .* (alpha .* t4 + t1 .* t4 + alpha) .* t1 .* F ./ 0.2e1;
I = t16(:,:,2) - t16(:,:,1);

end
