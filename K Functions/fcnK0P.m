function I = fcnK0P(S, T, u, alpha, F, tol)
t1 = F .^ 2;t4 = sqrt(t1);t9 = -0.1e1 ./ t1 .^ 2 ./ F ./ t4 ./ 0.6e1;
I = t9(:,:,2) - t9(:,:,1);

end
