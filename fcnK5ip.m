function I = fcnK5ip(S, t, u, alpha, F, tol)

F(abs(F(:,:,1)) < tol,:,1) = sign(F(abs(F(:,:,1)) < tol,:,1)).*tol;
F(abs(F(:,:,2)) < tol,:,2) = sign(F(abs(F(:,:,2)) < tol,:,2)).*tol;
t1 = F .^ 2;
t3 = F .* t;
t5 = sqrt(t1 .* S + t3 + u);
t12 = t .^ 2;
t17 = -0.2e1 ./ (0.4e1 .* S .* u - t12) .* (t3 + 0.2e1 .* u) ./ t5;

I = t17(:,:,2) - t17(:,:,1);

end
