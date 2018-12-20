function I = fcnK7ip(S, t, u, alpha, F, tol)

F(abs(F(:,:,1)) < tol,:,1) = sign(F(abs(F(:,:,1)) < tol,:,1)).*tol;
F(abs(F(:,:,2)) < tol,:,2) = sign(F(abs(F(:,:,2)) < tol,:,2)).*tol;
t1 = F .^ 2;
t2 = (t .^ 2);
t8 = u .^ 2;
t11 = sqrt(S);
t14 = S .^ 2;
t15 = t11 .* t14;
t19 = F .* t;
t24 = S .* u;
t27 = log(0.2e1);
t32 = sqrt(t1 .* S + t19 + u);
t36 = log(0.2e1 .* S .* F + 0.2e1 .* t11 .* t32 + t);
t37 = log(S);
t55 = -0.24e2 ./ (0.16e2 .* t24 - (4 .* t2)) ./ t32 ./ t15 .* (t11 .* S .* (t2 .* t1 ./ 0.6e1 - 0.5e1 ./ 0.3e1 .* u .* t .* F - 0.4e1 ./ 0.3e1 .* t8) - 0.2e1 ./ 0.3e1 .* u .* t15 .* t1 + t .* (t11 .* (t19 + u) .* t ./ 0.2e1 + t32 .* (-t27 + t36 - t37 ./ 0.2e1) .* (t24 - t2 ./ 0.4e1)));

I = t55(:,:,2) - t55(:,:,1);

end
