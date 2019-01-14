function I = fcnK1H(S, T, u, alpha, F, tol)
t1 = F .^ 2;t3 = S .* t1;t4 = sqrt(t3);t9 = -0.1e1 ./ t1 ./ t4 ./ t3 ./ 0.5e1;
I = t9(:,:,2) - t9(:,:,1);

end
