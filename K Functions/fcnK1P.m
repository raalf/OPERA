function I = fcnK1P(S, T, u, alpha, F, tol)
t1 = F .^ 2;t3 = sqrt(t1);t8 = -0.1e1 ./ t1 .^ 2 ./ t3 ./ 0.5e1;
I = t8(:,:,2) - t8(:,:,1);

end
