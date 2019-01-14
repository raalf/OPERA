function I = fcnK4N(S, T, u, alpha, F, tol)
t1 = F .^ 2;t3 = sqrt(t1 + u);t7 = F ./ u ./ t3;
I = t7(:,:,2) - t7(:,:,1);

end
