function I = fcnK5B(S, t, u, alpha, F, tol)
t1 = F .^ 2;t3 = t .* F;t5 = sqrt(t1 .* S + t3 + u);t12 = t .^ 2;t17 = -0.2e1 ./ (0.4e1 .* S .* u - t12) .* (t3 + 0.2e1 .* u) ./ t5;
I = t17(:,:,2) - t17(:,:,1);

end
