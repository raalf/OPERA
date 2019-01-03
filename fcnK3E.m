function I = fcnK3E(S, t, u, alpha, F, tol)
t1 = (F .^ 2);t5 = t1 + u;t6 = t5 .^ 2;t7 = sqrt(t5);t12 = -(5 .* t1 + 2 .* u) ./ t7 ./ t6 ./ 0.15e2;
I = t12(:,:,2) - t12(:,:,1);

end
