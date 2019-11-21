function [int_circ] = fcnINTCIRC2(matPLEX, matCOEFF, matDVEGRID)

xi_1 = permute(matPLEX(1,1,:),[3 2 1]);
xi_2 = permute(matPLEX(2,1,:),[3 2 1]);
xi_3 = permute(matPLEX(3,1,:),[3 2 1]);

eta_1 = permute(matPLEX(1,2,:),[3 2 1]);
eta_2 = permute(matPLEX(2,2,:),[3 2 1]);
eta_3 = permute(matPLEX(3,2,:),[3 2 1]);

idx_flp = xi_3 < xi_1;

C = (eta_3 - eta_2)./(xi_3 - xi_2);
D_LE = eta_2 - ((xi_2.*(eta_3 - eta_2))./(xi_3 - xi_2));
E = (eta_3 - eta_1)./(xi_3 - xi_1);
D_TE = eta_1 - ((xi_1.*(eta_3 - eta_1))./(xi_3 - xi_1));

A_1 = matCOEFF(:,1); A_2 = matCOEFF(:,2); B_1 = matCOEFF(:,3);
B_2 = matCOEFF(:,4); C_2 = matCOEFF(:,5); C_3 = matCOEFF(:,6);

tmp = fcnINTCIRC(xi_1, xi_3, C, D_LE, E, D_TE, A_1, A_2, B_1, B_2, C_2, C_3);
tmp(~idx_flp) = -tmp(~idx_flp);

int_circ = sum(sum(tmp(matDVEGRID),1), 'all');

end