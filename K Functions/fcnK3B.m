function I = fcnK3B(S, T, u, alpha, F, tol)
t1 = 0.1e1 ./ u;t2 = F .^ 2;t4 = F .* T;t6 = sqrt(t2 .* S + t4 + u);t7 = 0.1e1 ./ t6;t15 = T .^ 2;t21 = sqrt(u);t30 = log(0.1e1 ./ F .* (0.2e1 .* t6 .* t21 + t4 + 0.2e1 .* u));t32 = t7 .* t1 - t7 ./ (0.4e1 .* S .* u - t15) .* (0.2e1 .* S .* F + T) .* t1 .* T - t30 ./ t21 ./ u;
I = t32(:,:,2) - t32(:,:,1);

end
