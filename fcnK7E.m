function I = fcnK7E(S, t, u, alpha, F, tol)
t1 = (F .^ 2);t2 = (t1 .^ 2);t7 = (u .^ 2);t13 = t1 + u;t14 = t13 .^ 2;t15 = sqrt(t13);t19 = (5 .* t2 .* t1 + 40 .* t1 .* t7 + 30 .* t2 .* u + 16 .* t7 .* u) ./ t15 ./ t14 ./ 0.5e1;
I = t19(:,:,2) - t19(:,:,1);

end
