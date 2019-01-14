function I = fcnK6L(S, T, u, alpha, F, tol)
t1 = F .^ 2;t4 = sqrt(F .* T + t1);t10 = log(T ./ 0.2e1 + F + t4);t11 = -0.2e1 ./ t4 .* F + t10;
I = t11(:,:,2) - t11(:,:,1);

end
