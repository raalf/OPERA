function I = fcnK2ip(S, t, u, alpha, F, tol)

F(abs(F(:,:,1)) < tol,:,1) = sign(F(abs(F(:,:,1)) < tol,:,1)).*tol;
F(abs(F(:,:,2)) < tol,:,2) = sign(F(abs(F(:,:,2)) < tol,:,2)).*tol;
t1 = F .^ 2;
t3 = F .* t;
t5 = sqrt(t1 .* S + t3 + u);
t7 = S .^ 2;
t10 = S .* F;
t13 = t .^ 2;
t16 = sqrt(u);
t19 = u .^ 2;
t20 = t16 .* t19;
t29 = S .* u - t13 ./ 0.4e1;
t34 = log(0.2e1 .* t5 .* t16 + t3 + 0.2e1 .* u);
t35 = log(F);
t50 = 0.3e1 ./ 0.2e1 ./ F ./ t29 ./ t20 .* (t16 .* u .* (-0.4e1 ./ 0.3e1 .* t7 .* t1 - 0.5e1 ./ 0.3e1 .* t .* t10 + t13 ./ 0.6e1) - 0.2e1 ./ 0.3e1 .* t20 .* S + F .* (t16 .* (t10 + t) .* t ./ 0.2e1 + t5 .* (t34 - t35) .* t29) .* t) ./ t5;

I = t50(:,:,2) - t50(:,:,1);

end
