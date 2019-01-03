function I = fcnK4E(S, t, u, alpha, F, tol)
t1 = F .^ 2;t2 = t1 + u;t3 = t2 .^ 2;t4 = sqrt(t2);t9 = t1 .^ 2;t12 = 0.1e1 ./ u ./ t4 ./ t3 .* t9 .* F ./ 0.5e1;
I = t12(:,:,2) - t12(:,:,1);

end
