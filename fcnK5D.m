function I = fcnK5D(S, t, u, alpha, F, tol)
t1 = F .^ 2;t4 = sqrt(S .* t1 + u);t8 = -0.1e1 ./ S ./ t4;
I = t8(:,:,2) - t8(:,:,1);

end
