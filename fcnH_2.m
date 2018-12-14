function I = fcnH_2(C, B, A, X1, X2, tol)
% Solves an integral of the form,
%
%  X2
%   /
%  |            1
%  |  -------------------- dX
% /           2
%  X1 sqrt(C X  + B X + A)
%
% Using the method in Table of Integrals, Series, and Products by
% Gradshteyn & Ryzhik
X1(abs(X1) < tol) = sign(X2(abs(X1) < tol)).*tol;
X2(abs(X2) < tol) = sign(X1(abs(X2) < tol)).*tol;

len = size(C,1);

delta = 4.*A.*C - B.^2;

X = [X1 X2];

I1 = nan(len,2);

%% C > 0
idx = C > 0;
if any(idx)
    I1(idx,:) = (1./sqrt(C(idx))).*log( (2.*sqrt(C(idx).*sqrt(C(idx).*X(idx,:).^2 + B(idx).*X(idx,:) + A(idx)))+ 2.*C(idx).*X(idx,:) + B(idx))./(sqrt(delta(idx))));
end

%% C > 0 & delta > 0 & abs(B) > 0
idx = C > 0 & delta > 0 & abs(B) > tol;
if any(idx)
    I1(idx,:) = (1./sqrt(C(idx))).*asinh((2.*C(idx).*X(idx,:) + B(idx))./(sqrt(delta(idx))));
end
%% C < 0 & delta < 0 & abs(B) > 0
idx = C < 0 & delta < 0 & abs(B) > tol;
if any(idx)
    I1(idx,:) = (-1./sqrt(-C(idx))).*asin((2.*C(idx).*X(idx,:) + B(idx))./sqrt(-delta(idx)));
end

%% C > 0 & delta = 0 & abs(B) > 0
idx = C > 0 & abs(delta) < tol & abs(B) > tol;
if any(idx)
    I1(idx,:) = (1./sqrt(C(idx))).*log(2.*C(idx).*X(idx,:) + B(idx));
end

%% C > 0 & B == 0
idx = C > 0 & abs(B) < tol;
if any(idx)
    I1(idx,:) = (1./sqrt(C(idx))).*log(X(idx,:).*sqrt(C(idx)) + sqrt(A(idx) + C(idx).*X(idx,:).^2));
end

%% C > 0 & A > 0 & B == 0
idx = C < 0 & A > 0 & abs(B) < tol;
if any(idx)
    I1(idx,:) = (1./sqrt(-C(idx))).*asin(X(idx,:)).*sqrt(-C(idx)./A(idx)); % Is this right? Or is the last term in the arcsin?
end

%% C ~= 0, A == 0 & B == 0
idx = abs(C) > tol & abs(A) <= tol & abs(B) <= tol;
if any(idx)
   I1(idx,:) = (X(idx,:).*log(X(idx,:)))./(sqrt(C(idx,:).*X(idx,:).^2)); 
end

%%
I = I1(:,2) - I1(:,1);
I = real(I);
end