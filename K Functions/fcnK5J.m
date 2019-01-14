function I = fcnK5J(S, T, u, alpha, F, tol)
t1 = F .^ 2;t2 = F .* T;t4 = sqrt(t1 + t2 + u);t9 = T .^ 2;t14 = 0.2e1 ./ (t9 - 0.4e1 .* u) .* (t2 + 0.2e1 .* u) ./ t4;
I = t14(:,:,2) - t14(:,:,1);

end
