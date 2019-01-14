function I = fcnK4P(S, T, u, alpha, F, tol)
t1 = F .^ 2;t2 = sqrt(t1);t7 = -F ./ t2 ./ t1 ./ 0.2e1;
I = t7(:,:,2) - t7(:,:,1);

end
