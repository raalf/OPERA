function I = fcnK5D(S, T, u, alpha, F, tol)
t1 = F .^ 2;t8 = F .* T + S .* t1;t9 = sqrt(t8);t14 = 0.2e1 ./ t9 ./ t8 ./ T .* (F .* S + T) .* t1;
I = t14(:,:,2) - t14(:,:,1);

end
