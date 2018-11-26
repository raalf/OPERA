function I = fcnH_4(p, q, t1, t2)
% Solves an integral of the form, 
%
%  t2
%   /
%  |            t
%  |  --------------------- dt
% /     2            2
%  t1 (t  + p) sqrt(t  + q)
%
% Using the method in Table of Integrals, Series, and Products by
% Gradshteyn & Ryzhik (7th Edition) 
%
% Use u = sqrt(t.^2 + q) to transform into sign(u)./(u^2 + p - q)

u(:,1) = sqrt(t1.^2 + q);
u(:,2) = sqrt(t2.^2 + q);
a = (p - q);

I1 = fcnH_5(a, u);

% I = sign(u(:,1)).*I1(:,2) - sign(u(:,2)).*I1(:,1);
I = I1(:,2) - I1(:,1);
end