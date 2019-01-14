function I = fcnK6F(S, T, u, alpha, F, tol)
t3 = F .^ 2;t6 = sqrt(S .* t3 + u);t9 = sqrt(S);t14 = log(F .* t9 + t6);t16 = -0.1e1 ./ t6 ./ S .* F + t14 ./ t9 ./ S;
I = t16(:,:,2) - t16(:,:,1);

end
