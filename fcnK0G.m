function I = fcnK0G(S, t, u, alpha, F, tol)
t1 = F .^ 2;t4 = S .* t1;t5 = sqrt(t4);t10 = -0.1e1 ./ t1 ./ F ./ t5 ./ t4 ./ 0.6e1;
I = t10(:,:,2) - t10(:,:,1);

end
