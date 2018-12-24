function I = fcnH_7(M, N, alpha, F1, F2, tol)
% Solves an integral of the form,
%
%  F2
%   /
%  |                N + M F
%  |  arcsinh(------------------) dF
% /                  2
%  F1          sqrt(F  + alpha) 
%

%%
I = nan(size(M));
F1(abs(F1) < tol) = sign(F2(abs(F1) < tol)).*tol;
F2(abs(F2) < tol) = sign(F1(abs(F2) < tol)).*tol;

idx_in_plane = sqrt(alpha) < tol;
I(~idx_in_plane) = fcnH_7op(M(~idx_in_plane), N(~idx_in_plane), alpha(~idx_in_plane), F1(~idx_in_plane), F2(~idx_in_plane));
I(idx_in_plane) = fcnH_7ip(M(idx_in_plane), N(idx_in_plane), F1(idx_in_plane), F2(idx_in_plane));

I(abs(N) <= tol & abs(M) <= tol) = 0;
idx = abs(N) <= tol & abs(M) > tol & idx_in_plane;
I(idx) = F2(idx).*asinh(sign(F2(idx)).*M(idx)) - F1(idx).*asinh(sign(F1(idx)).*M(idx));

% I(abs(M) <= tol & abs(N) > tol) = fcnH_7C(N(abs(M) <= tol & abs(N) > tol), alpha(abs(M) <= tol & abs(N) > tol), F1(abs(M) <= tol & abs(N) > tol), F2(abs(M) <= tol & abs(N) > tol));
end

function I = fcnH_7C(N, alpha, F1, F2)
t1 = [F1 F2] .^ 2;
t2 = t1 + alpha;
t3 = sqrt(t2);
t6 = asinh(0.1e1 ./ t3 .* N);
t8 = N .^ 2;
t9 = t1 + t8 + alpha;
t12 = sqrt(0.1e1 ./ t2 .* t9);
t14 = sqrt(t9);
t16 = log([F1 F2] + t14);
t17 = sqrt(-alpha);
t19 = sign(N);
t20 = N .* t19;
t23 = t14 .* t20;
t24 = [F1 F2] .* t17;
t30 = log(-0.2e1 ./ (-[F1 F2] + t17) .* (t8 + t23 + alpha + t24));
t37 = log(0.2e1 ./ ([F1 F2] + t17) .* (t8 + t23 + alpha - t24));
t47 = t6 .* [F1 F2] + t19 ./ t17 ./ t14 .* (0.2e1 .* t20 .* t17 .* t16 + t30 .* alpha - t37 .* alpha) .* t3 .* t12 ./ 0.2e1;

if size(t47,2) == 2
    I = t47(:,2) - t47(:,1);   
else
    I = [];
end

end

function I = fcnH_7op(M, N, alpha, F1, F2)
t1 = M .* [F1 F2];
t3 = [F1 F2] .^ 2;
t4 = t3 + alpha;
t5 = sqrt(t4);
t8 = asinh(0.1e1 ./ t5 .* (t1 + N));
t10 = M .^ 2;
t11 = alpha .* t10;
t12 = M .* N;
t13 = sqrt(-alpha);
t15 = 0.2e1 .* t13 .* t12;
t16 = N .^ 2;
t18 = sqrt(-t11 + t15 + t16);
t21 = sqrt(-t11 - t15 + t16);
t24 = N .* t1;
t26 = t10 .* t3 + alpha + t16 + 0.2e1 .* t24 + t3;
t27 = sqrt(t26);
t29 = sqrt(t10 + 0.1e1);
t30 = t29 .* t27;
t36 = sqrt(0.1e1 ./ t4 .* t26);
t38 = [F1 F2] .* t10;
t40 = 0.1e1 ./ t29;
t42 = log(t40 .* (t38 + t12 + t30 + [F1 F2]));
t48 = t38 + t12 + [F1 F2];
t54 = log(0.1e1 ./ ([F1 F2] - t13) .* (t13 .* t48 + t27 .* t18 + alpha + t16 + t24));
t55 = log(0.2e1);
t59 = -t13 .* alpha .* M;
t60 = N .* alpha;
t71 = log(0.1e1 ./ ([F1 F2] + t13) .* (-t13 .* t48 + t27 .* t21 + alpha + t16 + t24));
t90 = -0.1e1 ./ t13 .* t40 ./ t21 ./ t18 ./ t27 .* (-0.2e1 .* t30 .* t21 .* t13 .* t18 .* t8 .* [F1 F2] + (t21 .* (-0.2e1 .* t13 .* N .* t18 .* t42 + (t59 - t60) .* (t54 + t55) .* t29) + t18 .* (t59 + t60) .* t29 .* (t71 + t55)) .* t36 .* t5) ./ 0.2e1;
if size(t90,2) == 2
    I = t90(:,2) - t90(:,1);   
else
    I = [];
end
end

function I = fcnH_7ip(M, N, F1, F2)
t1 = [F1 F2] .^ 2;
t2 = [F1 F2] .* M;
t4 = sign([F1 F2]);
t6 = 0.1e1 ./ [F1 F2];
t8 = asinh(t6 .* t4 .* (t2 + N));
t10 = M .^ 2;
t14 = N .^ 2;
t15 = 0.2e1 .* N .* t2 + t10 .* t1 + t1 + t14;
t18 = sqrt(0.1e1 ./ t1 .* t15);
t20 = sqrt(t10 + 0.1e1);
t24 = sqrt(t15);
t28 = 0.1e1 ./ t20;
t30 = log(t28 .* ([F1 F2] .* t10 + N .* M + t24 .* t20 + [F1 F2]));
t38 = t6 ./ t18 .* t28 .* (t24 .* t4 .* t30 .* N + t20 .* t18 .* t8 .* t1);
if size(t38,2) == 2
    I = t38(:,2) - t38(:,1);   
else
    I = [];
end
end
