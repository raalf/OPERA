function I = fcnK5H(S, T, u, alpha, F, tol)
t1 = F .^ 2;t2 = S .* t1;t3 = sqrt(t2);t7 = -t1 ./ t3 ./ t2;
I = t7(:,:,2) - t7(:,:,1);

end
