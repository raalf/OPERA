function I = fcnK4ip(S, t, u, alpha, F, tol)

F(abs(F(:,:,1)) < tol,:,1) = sign(F(abs(F(:,:,1)) < tol,:,1)).*tol;
F(abs(F(:,:,2)) < tol,:,2) = sign(F(abs(F(:,:,2)) < tol,:,2)).*tol;
t5 = F .^ 2;
t9 = sqrt(F .* t + t5 .* S + u);
t14 = t .^ 2;
t17 = 0.1e1 ./ (0.4e1 .* S .* u - t14) ./ t9 .* (0.4e1 .* F .* S + 0.2e1 .* t);

I = t17(:,:,2) - t17(:,:,1);

end
