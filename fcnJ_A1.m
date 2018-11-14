function J_A1 = fcnJ_A1(x_m, y_m, z_m, xi_1, xi_3, C, D_LE)

t2 = (C .* x_m + D_LE - y_m);
t4 = abs(z_m);
t5 = t4 .* C;
t13 = ((x_m + z_m) .* C - y_m + D_LE) .* ((x_m - z_m) .* C - y_m + D_LE);
t15 = sqrt(-2.*i .* t5 .* t2 + t13);
t21 = sqrt(2.*i .* t5 .* t2 + t13);
t22 = -i .* t21;
t23 = C .^ 2;
t24 = (xi_1 .^ 2);
t26 = y_m - D_LE;
t30 = x_m .^ 2;
t33 = y_m .^ 2;
t35 = 2 .* D_LE .* y_m;
t36 = z_m .^ 2;
t37 = D_LE .^ 2;
t39 = sqrt((-2 .* C .* t26 .* xi_1 + t24 .* t23 - 2 .* x_m .* xi_1 + t24 + t30 + t33 - t35 + t36 + t37));
t41 = t23 .* xi_1;
t43 = -C .* t26;
t47 = -i .* t41 .* x_m;
t50 = C .* (x_m + xi_1) .* t26;
t51 = -i .* t33;
t53 = 2.*i .* D_LE .* y_m;
t54 = -i .* t36;
t55 = -i .* t37;
t57 = -i .* xi_1;
t61 = log(0.1e1 ./ (i .* x_m + t57 - t4) .* (t39 .* t22 + (t4 .* (t41 + t43 - x_m + xi_1)) + t47 + i .* t50 + t51 + t53 + t54 + t55));
t62 = (xi_3 .^ 2);
t70 = sqrt((-2 .* C .* t26 .* xi_3 + t62 .* t23 - 2 .* x_m .* xi_3 + t30 + t33 - t35 + t36 + t37 + t62));
t72 = t23 .* xi_3;
t75 = -i .* xi_3;
t77 = t23 .* x_m .* t75;
t80 = C .* (x_m + xi_3) .* t26;
t85 = log(0.1e1 ./ (i .* x_m + t75 - t4) .* (t70 .* t22 + (t4 .* (t72 + t43 - x_m + xi_3)) + t77 + i .* t80 + t51 + t53 + t54 + t55));
t90 = -i .* t15;
t92 = C .* t26;
t99 = log(0.1e1 ./ (t4 + i .* x_m + t57) .* (t39 .* t90 + (t4 .* (-t41 + t92 + x_m - xi_1)) + t47 + i .* t50 + t51 + t53 + t54 + t55));
t107 = log(0.1e1 ./ (t4 + i .* x_m + t75) .* (t70 .* t90 + (t4 .* (-t72 + t92 + x_m - xi_3)) + t77 + i .* t80 + t51 + t53 + t54 + t55));
J_A1 = -0.1e1 ./ 0.2e1.*i ./ t4 ./ t21 .* (t15 .* (i .* t4 + x_m) .* (t61 - t85) + (i .* t4 - x_m) .* t21 .* (t99 - t107)) ./ t15;


end

