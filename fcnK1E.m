function I = fcnK1E(S, t, u, alpha, F, tol)
t1 = F .^ 2;t2 = t1 + u;t3 = t2 .^ 2;t4 = sqrt(t2);t8 = -0.1e1 ./ t4 ./ t3 ./ 0.5e1;
I = t8(:,:,2) - t8(:,:,1);

end
