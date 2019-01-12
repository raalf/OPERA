function I = fcnK2G(S, t, u, alpha, F, tol)
t2 = F .^ 2;t3 = S .* t2;t4 = sqrt(t3);t9 = -0.1e1 ./ F ./ t4 ./ t3 ./ 0.4e1;
I = t9(:,:,2) - t9(:,:,1);

end
