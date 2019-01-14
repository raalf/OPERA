function I = fcnK2T(S, T, u, alpha, F, tol)
t1 = (F .^ 2);t7 = t1 + u;t8 = t7 .^ 2;t9 = sqrt(t7);t12 = u .^ 2;t16 = t1 .* F .* (2 .* t1 + 5 .* u) ./ t9 ./ t8 ./ t12 ./ 0.15e2;
I = t16(:,:,2) - t16(:,:,1);

end
