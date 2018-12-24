function I = fcnK3D(S, t, u, alpha, F, tol)
t2 = F .^ 2;t5 = sqrt(S .* t2 + u);t8 = sqrt(u);t15 = log(0.2e1 ./ F .* (t5 .* t8 + u));t17 = 0.1e1 ./ t5 ./ u - t15 ./ t8 ./ u;
I = t17(:,:,2) - t17(:,:,1);

end
