function I = fcnK4G(S, T, u, alpha, F, tol)
t1 = F .^ 2;t6 = t1 .* S;t7 = sqrt(t6);t12 = -0.1e1 ./ t7 ./ t6 .* t1 .* F ./ (t1 + alpha) ./ 0.2e1;
I = t12(:,:,2) - t12(:,:,1);

end
