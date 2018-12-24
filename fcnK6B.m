function I = fcnK6B(S, t, u, alpha, F, tol)
t1 = 0.1e1 ./ S;t3 = F .^ 2;t7 = sqrt(t .* F + t3 .* S + u);t8 = 0.1e1 ./ t7;t10 = S .^ 2;t11 = 0.1e1 ./ t10;t15 = t .^ 2;t21 = t8 ./ (0.4e1 .* S .* u - t15);t28 = sqrt(S);t37 = log(0.1e1 ./ t28 .* (t ./ 0.2e1 + S .* F) + t7);t39 = -t8 .* t1 .* F + t8 .* t11 .* t ./ 0.2e1 + F .* t21 .* t1 .* t15 + t21 .* t11 .* t15 .* t ./ 0.2e1 + t37 ./ t28 ./ S;
I = t39(:,:,2) - t39(:,:,1);

end
