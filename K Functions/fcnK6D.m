function I = fcnK6D(S, T, u, alpha, F, tol)
t3 = F .^ 2;t7 = sqrt(F .* T + S .* t3);t11 = sqrt(S);t20 = log(0.1e1 ./ t11 .* (T ./ 0.2e1 + F .* S) + t7);t22 = -0.2e1 ./ t7 ./ S .* F + t20 ./ t11 ./ S;
I = t22(:,:,2) - t22(:,:,1);

end
