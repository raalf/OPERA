function I = fcnK2F(S, t, u, alpha, F, tol)
t2 = F .^ 2;t3 = sqrt(t2);t8 = -0.1e1 ./ F ./ t3 ./ t2 ./ 0.4e1;
I = t8(:,:,2) - t8(:,:,1);

end
