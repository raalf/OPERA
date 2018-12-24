function I = fcnK7D(S, t, u, alpha, F, tol)
t1 = F .^ 2;t2 = S .* t1;t6 = sqrt(t2 + u);t9 = S .^ 2;t11 = 0.1e1 ./ t9 ./ t6 .* (t2 + 0.2e1 .* u);
I = t11(:,:,2) - t11(:,:,1);

end
