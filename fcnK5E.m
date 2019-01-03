function I = fcnK5E(S, t, u, alpha, F, tol)
t1 = (F .^ 2);t2 = (t1 .^ 2);t6 = (u .^ 2);t9 = t1 + u;t10 = t9 .^ 2;t11 = sqrt(t9);t16 = -(20 .* t1 .* u + 15 .* t2 + 8 .* t6) ./ t11 ./ t10 ./ 0.15e2;
I = t16(:,:,2) - t16(:,:,1);

end
