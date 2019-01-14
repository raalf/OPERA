function I = fcnK6H(S, T, u, alpha, F, tol)
t1 = F .^ 2;t2 = S .* t1;t3 = sqrt(t2);t8 = log(F);t9 = 0.1e1 ./ t3 ./ t2 .* t1 .* F .* t8;
I = t9(:,:,2) - t9(:,:,1);

end
