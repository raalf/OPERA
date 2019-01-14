function I = fcnK4L(S, T, u, alpha, F, tol)
t5 = T .^ 2;t8 = F .^ 2;t10 = F .* T + t8;t11 = sqrt(t10);t17 = -0.2e1 ./ t11 ./ t10 ./ t5 .* (2 .* F + T) .* F .* (F + T);
I = t17(:,:,2) - t17(:,:,1);

end
