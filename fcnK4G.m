function I = fcnK4G(S, t, u, alpha, F, tol)
t1 = F .^ 2;t2 = S .* t1;t3 = sqrt(t2);t8 = -F ./ t3 ./ t2 ./ 0.2e1;
I = t8(:,:,2) - t8(:,:,1);

end
