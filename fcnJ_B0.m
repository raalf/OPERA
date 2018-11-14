function J_B0 = fcnJ_B0(x_m, y_m, z_m, xi_1, xi_3, E, D_TE)

t2 = (E .* x_m + D_TE - y_m);
t4 = abs(z_m);
t5 = t4 .* E;
t13 = ((x_m + z_m) .* E + D_TE - y_m) .* ((x_m - z_m) .* E + D_TE - y_m);
t15 = sqrt(-2.*i .* t5 .* t2 + t13);
t21 = sqrt(2.*i .* t5 .* t2 + t13);
t23 = -i .* t21;
t24 = E .^ 2;
t25 = (xi_1 .^ 2);
t27 = -D_TE + y_m;
t31 = x_m .^ 2;
t34 = y_m .^ 2;
t36 = 2 .* D_TE .* y_m;
t37 = z_m .^ 2;
t38 = D_TE .^ 2;
t40 = sqrt((-2 .* E .* t27 .* xi_1 + t25 .* t24 - 2 .* x_m .* xi_1 + t25 + t31 + t34 - t36 + t37 + t38));
t42 = t24 .* xi_1;
t44 = -E .* t27;
t47 = -i .* x_m;
t48 = t42 .* t47;
t51 = E .* (x_m + xi_1) .* t27;
t53 = -i .* (t38 - t36 + t34 + t37);
t55 = -i .* xi_1;
t59 = log(0.1e1 ./ (i .* x_m + t55 - t4) .* (t40 .* t23 + (t4 .* (t42 + t44 + xi_1 - x_m)) + t48 + i .* t51 + t53));
t61 = (xi_3 .^ 2);
t69 = sqrt((-2 .* E .* t27 .* xi_3 + t61 .* t24 - 2 .* x_m .* xi_3 + t31 + t34 - t36 + t37 + t38 + t61));
t71 = t24 .* xi_3;
t74 = t71 .* t47;
t77 = E .* (x_m + xi_3) .* t27;
t79 = -i .* xi_3;
t83 = log(0.1e1 ./ (i .* x_m + t79 - t4) .* (t69 .* t23 + (t4 .* (t71 + t44 + xi_3 - x_m)) + t74 + i .* t77 + t53));
t85 = -i .* t15;
t87 = E .* t27;
t90 = -i .* t38;
t92 = 2.*i .* D_TE .* y_m;
t93 = -i .* t34;
t94 = -i .* t37;
t99 = log(0.1e1 ./ (t4 + i .* x_m + t55) .* (t40 .* t85 + (t4 .* (-t42 + t87 - xi_1 + x_m)) + t48 + i .* t51 + t90 + t92 + t93 + t94));
t107 = log(0.1e1 ./ (t4 + i .* x_m + t79) .* (t69 .* t85 + (t4 .* (-t71 + t87 - xi_3 + x_m)) + t74 + i .* t77 + t90 + t92 + t93 + t94));
J_B0 = -0.1e1 ./ 0.2e1.*i ./ t4 .* (t15 .* t59 - t15 .* t83 - (t99 - t107) .* t21) ./ t21 ./ t15;

end

