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

X = [X1 X2];

I1 = nan(len,2);

%% C > 0 & delta > 0 & abs(B) > 0
idx = C > 0 & delta > 0 & abs(B) > zer;
if any(idx)
    I1(idx,:) = (1./sqrt(C(idx))).*asinh((2.*C(idx).*X(idx,:) + B(idx))./(sqrt(delta(idx))));
end
%% C < 0 & delta < 0 & abs(B) > 0
idx = C < 0 & delta < 0 & abs(B) > zer;
if any(idx)
    I1(idx,:) = (-1./sqrt(-C(idx))).*asin((2.*C(idx).*X(idx,:) + B(idx))./sqrt(-delta(idx)));
end

%% C > 0 & delta = 0 & abs(B) > 0
idx = C > 0 & abs(delta) < zer & abs(B) > zer;
if any(idx)
    I1(idx,:) = (1./sqrt(C(idx))).*log(2.*C(idx).*X(idx,:) + B(idx));
end

%% C > 0 & B == 0
idx = C > 0 & abs(B) < zer;
if any(idx)
    I1(idx,:) = (1./sqrt(C(idx))).*log(X(idx,:).*sqrt(C(idx)) + sqrt(A(idx) + C(idx).*X(idx,:).^2));
end

%% C > 0 & A > 0 & B == 0
idx = C < 0 & A > 0 & abs(B) < zer;
if any(idx)
    I1(idx,:) = (1./sqrt(-C(idx))).*asin(X(idx,:)).*sqrt(-C(idx)./A(idx)); % Is this right? Or is the last term in the arcsin?
end

%%
I = I1(:,2) - I1(:,1);
I = real(I);
end