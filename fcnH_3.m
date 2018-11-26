function I = fcnH_3(p, q, t1, t2)
% Solves an integral of the form, 
%
%  t2
%   /
%  |            1
%  |  --------------------- dt
% /     2            2
%  t1 (t  + p) sqrt(t  + q)
%
% Using the method in Table of Integrals, Series, and Products by
% Gradshteyn & Ryzhik (7th Edition) (2.172)
%
% Use 1/sqrt(t^2 + q) to transform into 1/(v^2 - p/(p - q))

v(:,1) = t1./(sqrt(t1.^2 + q));
v(:,2) = t2./(sqrt(t2.^2 + q));
a = -p./(p - q);

I1 = -fcnH_5(a, v);

I = I1(:,2) - I1(:,1);
end