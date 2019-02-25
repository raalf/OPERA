function I = fcnK6P(S, T, u, alpha, F, tol)
t1 = F .^ 2;
t2 = sqrt(t1);
t7 = log(F);
% t7 = sqrt(real(t7).^2 + imag(t7).^2);

t8 = 0.1e1 ./ t2 .* F .* t7;

I = t8(:,:,2) - t8(:,:,1);
I = sqrt(real(I).^2 + imag(I).^2);
end
