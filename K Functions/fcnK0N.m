function I = fcnK0N(S, T, u, alpha, F, tol)
t1 = F .^ 2;
t2 = (t1 .^ 2);
t6 = u .^ 2;
t9 = sqrt(t1 + u);
t18 = 0.1e1 ./ t6 ./ u ./ t1 ./ F ./ t9 .* (0.4e1 .* u .* t1 + (8 .* t2) - t6) ./ 0.3e1;

I = t18(:,:,2) - t18(:,:,1);

end
