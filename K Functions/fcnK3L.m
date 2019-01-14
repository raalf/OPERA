function I = fcnK3L(S, T, u, alpha, F, tol)
t2 = (F .^ 2);t4 = (F .* T);t6 = (T .^ 2);t11 = t2 + t4;t12 = sqrt(t11);t17 = 0.2e1 ./ 0.3e1 ./ t12 ./ t11 ./ t6 ./ T .* (8 .* t2 + 4 .* t4 - t6) .* (F + T);
I = t17(:,:,2) - t17(:,:,1);

end
