function I = fcnK7G(S, t, u, alpha, F, tol)
t1 = F .^ 2;t2 = S .* t1;t3 = sqrt(t2);t6 = t1 .^ 2;t7 = 0.1e1 ./ t3 ./ t2 .* t6;
I = t7(:,:,2) - t7(:,:,1);

end
