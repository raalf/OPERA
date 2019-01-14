function I = fcnK2G(S, T, u, alpha, F, tol)
t1 = F .^ 2;t3 = log(F);t6 = t1 + alpha;t7 = log(t6);t14 = t1 .* S;t15 = sqrt(t14);t20 = alpha .^ 2;t24 = 0.1e1 ./ t20 ./ t6 ./ t15 ./ t14 .* (0.2e1 .* alpha .* t3 - alpha .* t7 + 0.2e1 .* t1 .* t3 - t1 .* t7 + alpha) .* t1 .* F ./ 0.2e1;
I = t24(:,:,2) - t24(:,:,1);

end
