function I = fcnK4F(S, T, u, alpha, F, tol)
t1 = F .^ 2;t4 = sqrt(S .* t1 + u);t8 = F ./ u ./ t4;
I = t8(:,:,2) - t8(:,:,1);

end
