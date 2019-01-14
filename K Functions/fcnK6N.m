function I = fcnK6N(S, T, u, alpha, F, tol)
t1 = F .^ 2;t3 = sqrt(t1 + u);t7 = log(F + t3);t8 = -0.1e1 ./ t3 .* F + t7;
I = t8(:,:,2) - t8(:,:,1);

end
