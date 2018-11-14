function J_C = fcnJ_C(x_m, y_m, z_m, xi_1, xi_3, E, D_TE)
t1 = sign(z_m);
t2 = (E .^ 2);
t3 = 1 + t2;
t4 = t3 .^ (-0.1e1 ./ 0.2e1);
t5 = t3 .* t4;
t6 = (y_m - D_TE);
t7 = (D_TE .* y_m);
t8 = (t6 .* E);
t9 = (D_TE .^ 2);
t10 = (x_m .^ 2);
t11 = z_m .^ 2;
t12 = (y_m .^ 2);
t13 = t2 .* t10;
t14 = t11 .* t2;
t15 = 2.*i .* E .* t1;
t16 = 2;
t17 = t16 .* (t8 .* x_m + t7);
t18 = -t17 + t15 .* y_m .* z_m + t12 + t13 - t14 + t9 + -2.*i .* D_TE .* E .* t1 .* z_m + -2.*i .* t2 .* t1 .* x_m .* z_m;
t19 = t18 .^ (-0.1e1 ./ 0.2e1);
t18 = t18 .* t19;
t20 = (xi_1 .^ 2);
t21 = -(t16 .* ((x_m + t8) .* xi_1 + t7)) + (t20 .* t3) + t10 + t11 + t12 + t9;
t22 = sqrt(t21);
t23 = i + i .* t2;
t24 = i .* x_m;
t25 = t23 .* xi_1;
t26 = -i .* t1;
t27 = t26 .* z_m;
t28 = i .* t1;
t29 = t28 .* z_m;
t2 = t2 .* x_m;
t30 = t2 .* xi_1;
t31 = t7 .* t16;
t32 = xi_1 - x_m + t29;
t32 = 0.1e1 ./ t32;
t32 = log((-t31 + t18 .* t22 + (t1 .* (t24 - t25) + z_m) .* z_m + t12 + t9 + t30 + (-xi_1 - x_m + t29) .* E .* y_m + (xi_1 + x_m + t27) .* E .* D_TE) .* t32);
t13 = (-t17 + -2.*i .* E .* t1 .* y_m .* z_m + t15 .* D_TE .* z_m + t9 + t13 - t14 + 2.*i .* t2 .* t1 .* z_m + t12);
t14 = (t13 .^ (-0.1e1 ./ 0.2e1));
t13 = (t13 .* t14);
t15 = -i .* x_m;
t17 = ((xi_1 + x_m + t29) .* E);
t33 = xi_1 - x_m + t27;
t33 = 0.1e1 ./ t33;
t17 = log((-t31 + t13 .* t22 + (t1 .* (t25 + t15) + z_m) .* z_m + t12 + t30 + t9 - t17 .* y_m + t17 .* D_TE) .* t33);
t25 = log(0.2e1);
t30 = (E .* xi_1);
t22 = log((D_TE - y_m + t30) .* E + t22 .* t5 - x_m + xi_1);
t33 = (t16 .* x_m);
t34 = 0.1e1 ./ z_m;
t35 = log(t3);
t36 = D_TE - y_m;
t37 = -i .* t18;
t38 = i .* D_TE;
t39 = i .* t18;
t40 = (t32 + t25) .* t13;
t41 = t13 .* t18;
t42 = t24 .* E;
t43 = -i .* y_m;
t44 = 2.*i .* t18;
t45 = t32 + t25;
t46 = -t17 - t25;
t47 = t26 .* t13;
t48 = t5 .* z_m;
t49 = t48 .* x_m .* E;
t50 = (xi_3 .^ 2);
t3 = ((t16 .* ((-x_m - t8) .* xi_3 - t7)) + (t3 .* t50) + t10 + t11 + t12 + t9);
t7 = sqrt(t3);
t8 = (E .* xi_3);
t51 = log((D_TE - y_m + t8) .* E + (t5 .* t7) - x_m + xi_3);
t23 = t23 .* xi_3;
t2 = t2 .* xi_3;
t52 = xi_3 - x_m + t29;
t52 = 0.1e1 ./ t52;
t24 = log((t18 .* t7 + (t1 .* (t24 - t23) + z_m) .* z_m + t12 + t9 + t2 + (-xi_3 - x_m + t29) .* E .* y_m + (xi_3 + x_m + t27) .* E .* D_TE - t31) .* t52);
t29 = (xi_3 + x_m + t29) .* E;
t27 = xi_3 - x_m + t27;
t27 = 0.1e1 ./ t27;
t2 = log(((t13 .* t7) + (t1 .* (t23 + t15) + z_m) .* z_m + t12 + t2 + t9 - t29 .* y_m + t29 .* D_TE - t31) .* t27);
t7 = (t25 + t2);
t9 = (t6 .* z_m);
t12 = -i .* E .* t11;
t23 = (t6 .* t13);
t2 = (t16 .* t49 .* (t13 .* (t25 + t24) + t18 .* (-t25 - t2)) + t1 .* (t18 .* (t15 .* t13 .* (2 .* E .* t51 + t5 .* log(-t33 .* xi_3 + t10 + t11 + t50)) + t12 .* t5 .* t2 + 2.*i .* t48 .* t13 .* atan((xi_3 - x_m) .* t34)) + t13 .* t5 .* ((t42 + t43) .* x_m + t12) .* t24) + t5 .* (t18 .* (t1 .* (t38 .* t7 + t43 .* t7) .* x_m + t9 .* t2 + t26 .* E .* t11 .* t25 + t28 .* E .* t7 .* t10 + 2.*i .* t1 .* (log(sqrt(t3) - D_TE + y_m - t8) - 1) .* t13 .* xi_3) + t13 .* t25 .* t1 .* ((t42 + -i .* t6) .* x_m + t12) + t13 .* (t28 .* D_TE .* x_m - t9) .* t24) + t41 .* t1 .* (t42 + t38 + t43) .* t35 + t23 .* t44 .* t1 .* t51 + t48 .* (t18 .* t6 - t23) .* t25);
J_C = -0.1e1 ./ 0.2e1.*i .* t1 .* t19 .* t14 .* t4 .* (t16 .* t49 .* (t18 .* (-t17 - t25) + t40) + t1 .* (t5 .* ((t13 .* (i .* t36 .* t25 + t37 .* log(-t33 .* xi_1 + t10 + t11 + t20) + t38 .* t32) + t39 .* t25 .* t36 + t39 .* t36 .* t17) .* x_m + 2.*i .* xi_1 .* t18 .* t13 .* (log(-D_TE + y_m - t30 + sqrt(t21)) - 1) + 2.*i .* t41 .* z_m .* atan((xi_1 - x_m) .* t34) + E .* (t37 .* t17 + -i .* (t18 + t13) .* t25) .* t11 + E .* (i .* (t17 + t25) .* t18 + i .* t40) .* t10) + t44 .* t22 .* t13 .* (y_m - D_TE) + t41 .* (t42 + t38 + t43) .* t35) + t47 .* E .* t5 .* t32 .* t11 + t5 .* (t13 .* (D_TE .* t45 - t45 .* y_m) + t18 .* (D_TE .* t46 - t46 .* y_m)) .* z_m + t47 .* (2 .* E .* t18 .* t22 + y_m .* t5 .* t32) .* x_m + t2);
end