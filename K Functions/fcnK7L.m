function I = fcnK7L(S, T, u, alpha, F, tol)
t1 = F .^ 2;t4 = sqrt(F .* T + t1);t5 = 0.1e1 ./ t4;t12 = log(T ./ 0.2e1 + F + t4);t15 = t5 .* t1 + 0.3e1 .* F .* t5 .* T - 0.3e1 ./ 0.2e1 .* t12 .* T;
I = t15(:,:,2) - t15(:,:,1);

end
