function I = fcnK4B(S, T, u, alpha, F, tol)
t1 = F .^ 2;t5 = sqrt(F .* T + t1 .* S + u);t13 = T .^ 2;t17 = 0.2e1 ./ (0.4e1 .* S .* u - t13) .* (0.2e1 .* S .* F + T) ./ t5;
I = t17(:,:,2) - t17(:,:,1);

end
