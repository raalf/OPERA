function I = fcnK3P(S, T, u, alpha, F, tol)
t1 = F .^ 2;t2 = sqrt(t1);t6 = -0.1e1 ./ t2 ./ t1 ./ 0.3e1;
I = t6(:,:,2) - t6(:,:,1);

end
