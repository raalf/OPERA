function I = fcnK3D(S, T, u, alpha, F, tol)
t1 = (F .* S);t3 = (F .^ 2);t4 = (S .^ 2);t9 = (T .^ 2);t16 = F .* T + S .* t3;t17 = sqrt(t16);t22 = 0.2e1 ./ 0.3e1 ./ t17 ./ t16 ./ t9 ./ T .* (4 .* T .* t1 + 8 .* t4 .* t3 - t9) .* (t1 + T);
I = t22(:,:,2) - t22(:,:,1);

end
