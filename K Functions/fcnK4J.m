function I = fcnK4J(S, T, u, alpha, F, tol)
t1 = F .^ 2;t4 = sqrt(F .* T + t1 + u);t9 = T .^ 2;t15 = -0.2e1 ./ (t9 - 0.4e1 .* u) .* (0.2e1 .* F + T) ./ t4;
I = t15(:,:,2) - t15(:,:,1);

end
