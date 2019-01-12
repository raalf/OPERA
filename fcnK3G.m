function I = fcnK3G(S, t, u, alpha, F, tol)
t1 = F .^ 2;t2 = S .* t1;t3 = sqrt(t2);t7 = -0.1e1 ./ t3 ./ t2 ./ 0.3e1;
I = t7(:,:,2) - t7(:,:,1);

end
