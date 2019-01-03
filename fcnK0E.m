function I = fcnK0E(S, t, u, alpha, F, tol)
t1 = (F .^ 2);t2 = (t1 .^ 2);t6 = (u .^ 2);t10 = t1 + u;t11 = t10 .^ 2;t12 = sqrt(t10);t19 = F .* (20 .* t1 .* u + 8 .* t2 + 15 .* t6) ./ t12 ./ t11 ./ t6 ./ u ./ 0.15e2;
I = t19(:,:,2) - t19(:,:,1);

end
