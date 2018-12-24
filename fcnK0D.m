function I = fcnK0D(S, t, u, alpha, F, tol)
t1 = (S .^ 2);t2 = F .^ 2;t3 = (t2 .^ 2);t6 = S .* t2;t9 = u .^ 2;t12 = sqrt(t6 + u);t21 = 0.1e1 ./ t9 ./ u ./ t2 ./ F ./ t12 .* ((8 .* t3 .* t1) + 0.4e1 .* u .* t6 - t9) ./ 0.3e1;
I = t21(:,:,2) - t21(:,:,1);

end
