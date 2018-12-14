function I = fcnH_1(M, N, S, T, u, alpha, F1, F2)
% Solves an integral of the form,
%
%  F2
%   /
%  |               N + M F
%  |  --------------------------------- dF
% /     2                  2
%  F1 (F  + alpha) sqrt(S F  + T F + u)
%
% Using the method of David Yeh
% https://books.google.ca/books?id=F7jiBQAAQBAJ&printsec=frontcover#v=onepage&q&f=false

%%
tol = 1e-10;
F1(abs(F1) < tol) = sign(F2(abs(F1) < tol)).*tol;
F2(abs(F2) < tol) = sign(F1(abs(F2) < tol)).*tol;

idx_in_plane = sqrt(alpha) < tol;
I = nan(size(M));
I(~idx_in_plane,1) = fcnH_1op(M(~idx_in_plane), N(~idx_in_plane), S(~idx_in_plane), T(~idx_in_plane), u(~idx_in_plane), alpha(~idx_in_plane), F1(~idx_in_plane), F2(~idx_in_plane));
I(idx_in_plane,1) = fcnH_1ip(M(idx_in_plane), N(idx_in_plane), S(idx_in_plane), T(idx_in_plane), u(idx_in_plane), F1(idx_in_plane), F2(idx_in_plane));
end

function I = fcnH_1op(M, N, S, T, u, alpha, F1, F2)
% Out of plane
t1 = S .* alpha;
t2 = sqrt(-alpha);
t3 = t2 .* T;
t5 = sqrt(-t1 - t3 + u);
t6 = 0.1e1 ./ t5;
t7 = 0.2e1 .* t1;
t8 = 0.2e1 .* t3;
t9 = 0.2e1 .* u;
t11 = 0.2e1 .* t2 .* S;
t13 = [F1 F2] + t2;
t14 = t13 .* (T - t11);
t15 = t13 .^ 2;
t18 = sqrt(t15 .* S - t1 + t14 - t3 + u);
t24 = log(0.1e1 ./ t13 .* (0.2e1 .* t18 .* t5 + t14 - t7 - t8 + t9));
t27 = 0.1e1 ./ t2;
t32 = sqrt(-t1 + t3 + u);
t33 = 0.1e1 ./ t32;
t35 = [F1 F2] - t2;
t36 = t35 .* (T + t11);
t37 = t35 .^ 2;
t40 = sqrt(t37 .* S - t1 + t3 + t36 + u);
t46 = log(0.1e1 ./ t35 .* (0.2e1 .* t40 .* t32 + t36 - t7 + t8 + t9));
t53 = N .* t24 .* t6 .* t27 ./ 0.2e1 - N .* t46 .* t33 .* t27 ./ 0.2e1 - M .* t24 .* t6 ./ 0.2e1 - M .* t46 .* t33 / 0.2e1;
if size(t53,2) == 2
    I = t53(:,2) - t53(:,1);
else
    I = [];
end
end

function I = fcnH_1ip(M, N, S, T, u, F1, F2)
% In plane
t1 = [F1 F2] .^ 2;
t3 = [F1 F2]  .* T;
t5 = sqrt(S .* t1 + t3 + u);
t7 = sqrt(u);
t13 = log(0.2e1 .* t5 .* t7 + t3 + 0.2e1 .* u);
t14 = log([F1 F2]);
t28 = -0.1e1 ./ [F1 F2]  ./ t7 ./ u .* (t7 .* t5 .* N + (u .* M - T .* N ./ 0.2e1) .* (t13 - t14) .* [F1 F2] );
if size(t28,2) == 2
    I = t28(:,2) - t28(:,1);
else
    I = [];
end
end