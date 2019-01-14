function I = fcnK6T(S, T, u, alpha, F, tol)
t1 = F .^ 2;
t2 = t1 .^ 2;
t4 = t1 + u;
t5 = t4 .^ 2;
t6 = sqrt(t4);

tmp = F + t6;
tmp(abs(tmp) < tol) = tol;
t19 = log(tmp);


t20 = -t2 .* F ./ t6 ./ t5 ./ 0.5e1 - t1 .* F ./ t6 ./ t4 ./ 0.3e1 - F ./ t6 + t19;

I = t20(:,:,2) - t20(:,:,1);

end
