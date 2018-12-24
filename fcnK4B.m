function I = fcnK4B(S, t, u, alpha, F, tol)
t1 = F .^ 2;t5 = sqrt(t .* F + t1 .* S + u);t13 = t .^ 2;t17 = 0.2e1 ./ (0.4e1 .* S .* u - t13) .* (0.2e1 .* S .* F + t) ./ t5;
I = t17(:,:,2) - t17(:,:,1);

end
