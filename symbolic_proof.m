clc
clear

syms z b00 H00 b10 H10 b01 H01 H02 H20 a00 a10 a01 H11

syms c1x c1y c1z c2x c2y c2z c3x c3y c3z

% Local reference frame assumption
c1 = [c1x c1y 0];
c2 = [c2x c2y 0];
c3 = [0 0 c3z]; 

soln_1 = z*(b00*H00 + b10*H10 + b01*H01)*c1 - z*(a00*H00 + a10*H10 + a01*H01)*c2 ...
    + (b00*H10 + b10*H20 + b01*H11 - a00*H01 - a10*H11 - a01*H02)*c3;

disp('dS1_x:')
disp(soln_1(1))

disp('dS1_y:')
disp(soln_1(2))

disp('dS1_z:')
disp(soln_1(3))

%%
dO1 = [ z*H00*c1x z*H10*c1x z*H01*c1x -z*H00*c2x -z*H10*c2x -z*H01*c2x; ...
        z*H00*c1y z*H10*c1y z*H01*c1y -z*H00*c2y -z*H10*c2y -z*H01*c2y; ...
        H10*c3z H20*c3z H11*c3z -H01*c3z -H11*c3z -H02*c3z ];

x1 = [b00; b10; b01; a00; a10; a01];

soln_2 = dO1*x1

%%
syms B2 B1 C2 A2 C2 A1

soln_3 = subs(soln_2, b00, B2);
soln_3 = subs(soln_3, b10, 2*B1);
soln_3 = subs(soln_3, b01, C2);
soln_3 = subs(soln_3, a00, A2);
soln_3 = subs(soln_3, a10, C2);
soln_3 = subs(soln_3, a01, 2*A1);

x2 = [A1; A2; B1; B2; C2];
[dO2, ~] = equationsToMatrix(soln_3,x2)


