clc
clear

main;

% J1 = 1/q

% Off of element
% assume( abs(z_m) > 0 | ...
%     (z_m == 0 & ...
%         ((x_m < xi_1 | x_m > xi_3) | ...
%         (((xi_1 < x_m) & (x_m < xi_3)) & (eta > (C*x_m + D_LE) | eta < (E*x_m + D_TE)) ))));
    
%     assumeAlso(...
%     (z_m == 0 & (x_m < xi_1)) |...
%     (z_m == 0 & (x_m > xi_3)) |...
%     (z_m == 0 & (y_m > (C*x_m + D_LE))) |...
%     (z_m == 0 & (y_m < (E*x_m + D_TE))))

% assumeAlso(z_m == 0 & (y_m > (C*x_m + D_LE)))
assumeAlso(z_m == 0)
j1 = int(1/(q),eta,eta_LE,eta_TE);
tp = children(j1);
assumeAlso(x_m < xi_1)
j1 = int(tp(1,1), xi, xi_1, xi_3);
matlabFunction(j1,'file','j1_si_out')