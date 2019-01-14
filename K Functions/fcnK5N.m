function I = fcnK5N(S, T, u, alpha, F, tol)
t1 = F .^ 2;t3 = sqrt(t1 + u);t5 = -0.1e1 ./ t3;
I = t5(:,:,2) - t5(:,:,1);

end
