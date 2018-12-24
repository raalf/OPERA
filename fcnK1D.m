function I = fcnK1D(S, t, u, alpha, F, tol)
t2 = F .^ 2;t7 = sqrt(S .* t2 + u);t8 = 0.1e1 ./ t7;t11 = u .^ 2;t16 = sqrt(u);t24 = log(0.2e1 ./ F .* (t7 .* t16 + u));t27 = -t8 ./ t2 ./ u ./ 0.2e1 - 0.3e1 ./ 0.2e1 .* t8 ./ t11 .* S + 0.3e1 ./ 0.2e1 .* t24 ./ t16 ./ t11 .* S;
I = t27(:,:,2) - t27(:,:,1);

end
