function I = fcnK7N(S, T, u, alpha, F, tol)
t1 = F .^ 2;t3 = sqrt(t1 + u);t7 = (t1 + 0.2e1 .* u) ./ t3;
I = t7(:,:,2) - t7(:,:,1);

end
