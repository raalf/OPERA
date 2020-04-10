function t30 = fcnINTCIRC3(xi_1, xi_3, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3)

t1 = (E .^ 2);
t5 = t1 .* A_1 + 2 .* C_2 .* E + B_1;
t6 = (xi_1 .^ 2);
t11 = 3 .* E .* (A_1 .* D_TE + A_2);
t13 = 3 .* C_2 .* D_TE;
t14 = 3 .* B_2;
t17 = (xi_3 .^ 2);
t21 = D_TE .^ 2;
t30 = ((xi_3 - xi_1) .* (t6 .* t5 + xi_1 .* (xi_3 .* t5 + t11 + t13 + t14) + t17 .* t5 + xi_3 .* (t11 + t13 + t14) + 3 .* t21 .* A_1 + 6 .* A_2 .* D_TE + 6 .* C_3)) ./ 0.6e1;

end

