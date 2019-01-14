function I = fcnK4O(S, T, u, alpha, F, tol)
t1 = F .^ 2;t6 = sqrt(t1);t11 = -0.1e1 ./ t6 .* F ./ (t1 + alpha) ./ 0.2e1;
I = t11(:,:,2) - t11(:,:,1);

end
