clc
clear

syms z b00 H00 b10 H10 b01 H01 H02 H20 a00 a10 a01 H11 ...
    c1x c1y c1z c2x c2y c2z c3x c3y c3z Gamma B2 B1 C2 A2 C3 ...
    A1 eta xsi zeta u v z c1 c2 c3

c1 = [c1x c1y 0];
c2 = [c2x c2y 0];
c3 = [0 0 c3z]; 

Gamma = 0.5*A1*xsi^2 + A2*xsi + 0.5*B2*eta^2 + B2*eta + C2*eta*xsi + C3;
gamma = gradient(Gamma, [xsi, eta, zeta]);

% r = u*c1 + v*c2 - z*c3;
% cross(r,gamma)/(norm(r)^3)

syms rho h
r = sqrt(rho^2 + h^2)


syms m n phi dphi a l1 l2
Hmn_phi = symfun(int(int(((cos(phi)^m)*(sin(phi)^n)*(rho^(m+n+1))/r^3), rho, 0, a*sec(phi)), phi, a/l1, a/l2), [m, n]);


% syms l a
% Hmn = subs(Hmn_phi(1,1), cos(phi), a/(sqrt(l^2 + a^2)));
% Hmn = subs(Hmn, sin(phi), l/(sqrt(l^2 + a^2)));
% Hmn = subs(Hmn, sin(2*phi), 2*l/(sqrt(l^2 + a^2)));
% Hmn = subs(Hmn, tan(phi), l/a);
% Hmn = subs(Hmn, dphi, a/(sqrt(l^2 + a^2)));

Hmn = int(Hmn, l)



%%
% syms c1x c1y c1z c2x c2y c2z c3x c3y c3z
% 
% Local reference frame assumption
% c1 = [c1x c1y 0];
% c2 = [c2x c2y 0];
% c3 = [0 0 c3z]; 
% 
% soln_1 = z*(b00*H00 + b10*H10 + b01*H01)*c1 - z*(a00*H00 + a10*H10 + a01*H01)*c2 ...
%     + (b00*H10 + b10*H20 + b01*H11 - a00*H01 - a10*H11 - a01*H02)*c3;
% 
% disp('dS1_x:')
% disp(soln_1(1))
% 
% disp('dS1_y:')
% disp(soln_1(2))
% 
% disp('dS1_z:')
% disp(soln_1(3))
% 
% %%
% dO1 = [ z*H00*c1x z*H10*c1x z*H01*c1x -z*H00*c2x -z*H10*c2x -z*H01*c2x; ...
%         z*H00*c1y z*H10*c1y z*H01*c1y -z*H00*c2y -z*H10*c2y -z*H01*c2y; ...
%         H10*c3z H20*c3z H11*c3z -H01*c3z -H11*c3z -H02*c3z ];
% 
% x1 = [b00; b10; b01; a00; a10; a01];
% 
% soln_2 = dO1*x1
% 
% %%
% syms B2 B1 C2 A2 C2 A1
% 
% soln_3 = subs(soln_2, b00, B2);
% soln_3 = subs(soln_3, b10, 2*B1);
% soln_3 = subs(soln_3, b01, C2);
% soln_3 = subs(soln_3, a00, A2);
% soln_3 = subs(soln_3, a10, C2);
% soln_3 = subs(soln_3, a01, 2*A1);
% 
% x2 = [A1; A2; B1; B2; C2];
% [dO2, ~] = equationsToMatrix(soln_3,x2)


