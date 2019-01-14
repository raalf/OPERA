function I = fcnK2N(S, T, u, alpha, F, tol)
t1 = F .^ 2;t5 = sqrt(t1 + u);t9 = u .^ 2;t13 = -0.1e1 ./ t5 .* (0.2e1 .* t1 + u) ./ t9 ./ F;
I = t13(:,:,2) - t13(:,:,1);

end
