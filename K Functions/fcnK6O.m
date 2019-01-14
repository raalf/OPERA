function I = fcnK6O(S, T, u, alpha, F, tol)
t1 = F .^ 2;t3 = t1 + alpha;t4 = log(t3);t10 = sqrt(t1);t15 = 0.1e1 ./ t10 ./ t3 .* (alpha .* t4 + t1 .* t4 + alpha) .* F ./ 0.2e1;
I = t15(:,:,2) - t15(:,:,1);

end
