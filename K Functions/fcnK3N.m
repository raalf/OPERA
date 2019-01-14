function I = fcnK3N(S, T, u, alpha, F, tol)
t2 = F .^ 2;t4 = sqrt(t2 + u);t7 = sqrt(u);t14 = log(0.2e1 ./ F .* (t7 .* t4 + u));t16 = 0.1e1 ./ t4 ./ u - t14 ./ t7 ./ u;
I = t16(:,:,2) - t16(:,:,1);

end
