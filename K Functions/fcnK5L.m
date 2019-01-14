function I = fcnK5L(S, T, u, alpha, F, tol)
t1 = F .^ 2;t6 = F .* T + t1;t7 = sqrt(t6);t12 = 0.2e1 ./ t7 ./ t6 ./ T .* (F + T) .* t1;
I = t12(:,:,2) - t12(:,:,1);

end
