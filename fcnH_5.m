function I1 = fcnH_5(a, v)
% Solves an integral of the form, 
%
%  v2
%   /
%  |            1
%  |  --------------------- dv
% /           2
%  v1       v    +   a
%
% Using the method in Table of Integrals, Series, and Products by
% Gradshteyn & Ryzhik (7th Edition) (2.124)
zer = 1e-14;
I1 = v.*nan;

delta = 4.*a;

%% a > 0
idx = a > 0;
if any(idx)
    I1(idx,:) = (2./sqrt(delta(idx))).*atan2(2.*v(idx,:), sqrt(delta(idx)));  
end

%% a < 0
idx = a < 0;
if any(idx)
    I1(idx,:) = (-2./sqrt(-delta(idx))).*atanh(2.*v(idx,:)./sqrt(-delta(idx)));
end

%% a == 0
idx = abs(a) < zer;
if any(idx)
    I1(idx,:) = -1./v(idx,:);  
end

%%
% I = I1(:,2) - I1(:,1);

end