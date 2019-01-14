function I = fcnK1N(S, T, u, alpha, F, tol)
t2 = F .^ 2;t6 = sqrt(t2 + u);t7 = 0.1e1 ./ t6;t10 = u .^ 2;t14 = sqrt(u);t21 = log(0.2e1 ./ F .* (t6 .* t14 + u));t24 = -t7 ./ t2 ./ u ./ 0.2e1 - 0.3e1 ./ 0.2e1 .* t7 ./ t10 + 0.3e1 ./ 0.2e1 .* t21 ./ t14 ./ t10;
I = t24(:,:,2) - t24(:,:,1);

end
