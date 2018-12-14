function I = fcnH_6(C, B, A, X1, X2, tol)
% Solves an integral of the form,
%
%  X2
%   /
%  |            X
%  |  -------------------- dX
% /           2
%  X1 sqrt(C X  + B X + A)
%
% Using the method in Table of Integrals, Series, and Products by
% Gradshteyn & Ryzhik, 7th Edition (page 96, 2.261)
X1(abs(X1) < tol) = sign(X2(abs(X1) < tol)).*tol;
X2(abs(X2) < tol) = sign(X1(abs(X2) < tol)).*tol;

X = [X1 X2];

I1 = (sqrt(C.*X.^2 + B.*X + A)./C);
I = I1(:,2) - I1(:,1) - (B./(2.*C)).*fcnH_2(C, B, A, X1, X2, tol);

I = real(I);
end