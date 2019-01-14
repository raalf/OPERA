function I = fcnK6J(S, T, u, alpha, F, tol)
t1 = F .^ 2;t4 = sqrt(F .* T + t1 + u);t5 = 0.1e1 ./ t4;t6 = t5 .* F;t9 = T .^ 2;t12 = 0.1e1 ./ (-t9 + 0.4e1 .* u);t21 = log(T ./ 0.2e1 + F + t4);t22 = -t6 + t5 .* T ./ 0.2e1 + t6 .* t12 .* t9 + t5 .* t12 .* t9 .* T ./ 0.2e1 + t21;
I = t22(:,:,2) - t22(:,:,1);

end
