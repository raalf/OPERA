function I = fcnK6H(S, T, u, alpha, F, tol)
% t1 = F .^ 2;
% t2 = S .* t1;
% t3 = sqrt(t2);
% t8 = log(abs(F));
% % t8 = sqrt(real(t8).^2 + imag(t8).^2);
% % t8 = real(t8) + imag(t8)
% 
% t9 = 0.1e1 ./ t3 ./ t2 .* t1 .* F .* t8;
% 
% I = t9(:,:,2) - t9(:,:,1);
% % I = sqrt(real(I).^2 + imag(I).^2);
% % I = abs(real(I)) + abs(imag(I));

t9 = log(F) ./ (S.^(3/2));
I = t9(:,:,2) - t9(:,:,1);
% I = real(I) + abs(imag(I))
end
