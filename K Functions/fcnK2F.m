function I = fcnK2F(S, T, u, alpha, F, tol)
t1 = F .^ 2;t2 = S .* t1;t6 = sqrt(t2 + u);t10 = u .^ 2;t14 = -0.1e1 ./ t10 ./ F ./ t6 .* (0.2e1 .* t2 + u);
I = t14(:,:,2) - t14(:,:,1);

end
