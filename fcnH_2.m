function I = fcnH_2(C, B, A, X1, X2)
% Solves an integral of the form, 
%
%  X2
%   /
%  |            1
%  |  -------------------- dF
% /           2
%  X1 sqrt(C X  + B X + A)
%
% Using the method in Table of Integrals, Series, and Products by
% Gradshteyn & Ryzhik

zer = 1e-14;

len = size(C,1);

delta = 4.*A.*C - B.^2;

I1 = nan(len,1);
I2 = nan(len,1);

%% C > 0 & delta > 0 & abs(B) > 0
idx = C > 0 & delta > 0 & abs(B) > zer;
I1(idx) = (1./sqrt(C(idx))).*asinh((2.*C(idx).*X1(idx) + B(idx))./(sqrt(delta(idx))));
I2(idx) = (1./sqrt(C(idx))).*asinh((2.*C(idx).*X2(idx) + B(idx))./(sqrt(delta(idx))));

%% C < 0 & delta < 0 & abs(B) > 0
idx = C < 0 & delta < 0 & abs(B) > zer;
I1(idx) = (-1./sqrt(-C(idx))).*asin((2.*C(idx).*X1(idx) + B(idx))./sqrt(-delta(idx)));
I2(idx) = (-1./sqrt(-C(idx))).*asin((2.*C(idx).*X2(idx) + B(idx))./sqrt(-delta(idx)));

%% C > 0 & delta = 0 & abs(B) > 0
idx = C > 0 & abs(delta) < zer & abs(B) > zer;
I1(idx) = (1./sqrt(C(idx))).*log(2.*C(idx).*X1(idx) + B(idx)); 
I2(idx) = (1./sqrt(C(idx))).*log(2.*C(idx).*X2(idx) + B(idx)); 

%% C > 0 & B == 0
idx = C > 0 & abs(B) < zer;
I1(idx) = (1./sqrt(C(idx))).*log(X1(idx).*sqrt(C(idx)) + sqrt(A(idx) + C(idx).*X1(idx).^2));
I2(idx) = (1./sqrt(C(idx))).*log(X2(idx).*sqrt(C(idx)) + sqrt(A(idx) + C(idx).*X2(idx).^2));

%% C > 0 & A > 0 & B == 0
idx = C < 0 & A > 0 & abs(B) < zer;
I1(idx) = (1./sqrt(-C(idx))).*asin(X1(idx)).*sqrt(-C(idx)./A(idx)); % Is this right? Or is the last term in the arcsin?
I2(idx) = (1./sqrt(-C(idx))).*asin(X2(idx)).*sqrt(-C(idx)./A(idx)); % Is this right? Or is the last term in the arcsin?

%%
I = I2 - I1;
I = real(I);
end