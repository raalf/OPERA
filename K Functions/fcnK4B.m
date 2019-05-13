function I = fcnK4B(S, T, u, alpha, F, tol)
t5 = F .^ 2;t9 = sqrt(F .* T + t5 .* S + u);t14 = T .^ 2;t17 = 0.1e1 ./ (0.4e1 .* S .* u - t14) ./ t9 .* (0.4e1 .* S .* F + 0.2e1 .* T);
I = t17(:,:,2) - t17(:,:,1);

end
