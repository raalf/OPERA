function [u_out] = fcnKJU(xi_1, xi_3, eta_1, eta_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3, u, v, w, rho, lim)
u_1 = u(:,1); u_2 = u(:,2); u_3 = u(:,3);
v_1 = v(:,1); v_2 = v(:,2); v_3 = v(:,3);
w_1 = w(:,1); w_2 = w(:,2); w_3 = w(:,3);

t1 = (E .^ 2);
t6 = (t1 .* A_1 + 2 .* C_2 .* E + B_1) .* w_1;
t7 = xi_1 .^ 2;
t8 = (t7 .^ 2);
t14 = A_1 .* w_2 .* t1;
t16 = A_1 .* D_TE + A_2;
t21 = E .* (2 .* w_2 .* C_2 + 2 .* w_1 .* t16);
t24 = 2 .* w_1 .* (C_2 .* D_TE + B_2);
t25 = w_2 .* B_1;
t29 = (xi_3 .^ 2);
t32 = t14 + t21 + t24 + t25;
t36 = 0.4e1 ./ 0.3e1 .* A_1 .* w_3 .* t1;
t37 = D_TE .* w_2;
t39 = w_2 .* A_2;
t40 = w_3 .* C_2;
t42 = 0.8e1 ./ 0.3e1 .* E .* (A_1 .* t37 + t39 + t40);
t44 = D_TE .^ 2;
t45 = t44 .* A_1;
t47 = A_2 .* D_TE;
t50 = w_1 .* (0.8e1 ./ 0.3e1 .* C_3 + 0.4e1 ./ 0.3e1 .* t45 + 0.8e1 ./ 0.3e1 .* t47);
t52 = 0.4e1 ./ 0.3e1 .* w_3 .* B_1;
t54 = 0.8e1 ./ 0.3e1 .* w_2 .* B_2;
t56 = 0.8e1 ./ 0.3e1 .* C_2 .* t37;
t59 = (t29 .* xi_3);
t63 = t36 + t42 + t50 + t52 + t54 + t56;
t67 = 4 .* E .* t16 .* w_3;
t70 = 2 .* A_1 .* t44 .* w_2;
t72 = 4 .* D_TE .* (t39 + t40);
t74 = 0.4e1 .* w_2 .* C_3;
t76 = 4 .* w_3 .* B_2;
t79 = (t29 .^ 2);
t95 = rho .* (-eta_3 + eta_1) .* (0.4e1 ./ 0.5e1 .* t8 .* t6 + t7 .* xi_1 .* (0.4e1 ./ 0.5e1 .* xi_3 .* t6 + t14 + t21 + t24 + t25) + t7 .* (0.4e1 ./ 0.5e1 .* t29 .* t6 + (xi_3 .* t32) + t36 + t42 + t50 + t52 + t54 + t56) + xi_1 .* (0.4e1 ./ 0.5e1 .* t59 .* t6 + (t29 .* t32) + xi_3 .* t63 + t67 + t70 + t72 + t74 + t76) + 0.4e1 ./ 0.5e1 .* t79 .* t6 + (t59 .* t32) + t29 .* t63 + xi_3 .* (t67 + t70 + t72 + t74 + t76) + 0.4e1 .* (t45 + (2 .* t47) + 0.2e1 .* C_3) .* w_3) ./ 0.8e1;

u_out = t95;
end