clc
clear
%runs a single triangle wake element from freewake

% P - point that is being induced
% xo - reference point of inducing DVE
testDVE = [0.5 0.5 0; 0.5 0.5 0; -0.5 0.5 0; -0.5 -0.5 0];

nu = 0;
eps = 0;
psi = 0;
phiLE = 0;
phiTE = 45;

xo = [-0.25 0 0];


eta = 0.5;
xsi = 0.25;

DVE_type = 1; % wake DVE, no filament

Temp.DBL_EPS = 1e-14;
singfct = 0;

% [a3, b3, c3] = fcnDVEInduction(Temp, FP, xo, nu, eps, phiLE, phiTE, psi, eta, xsi, DVE_type, singfct);
% D = [c3' b3' a3'];

%%

granularity = 0.1;

%COEFF = [A1 A2 B1 B2 C]';
COEFF = [0 1 0 0 0]'; %only works with the first sheet 

% D = INFLUENCECOEFF2(testDVE, DVECT, FP);
count = 1;

for x_p = 0.4
    for y_p = -3:granularity:3
        for z_p = -3:granularity:3
            
            FP = [x_p y_p z_p];
            
            [~, cA2, cA1] = fcnDVEInduction(Temp, FP, xo, nu, eps, phiLE, phiTE, psi, eta, xsi, DVE_type, singfct);
            D = [cA1' cA2'];
%             A1 = c1
%             A2 = b1
%             A3 = a1
            [~, cB2, cB1] = fcnDVEInduction(Temp, FP, xo, nu, eps, phiLE, phiTE, psi+90, eta, xsi, DVE_type, singfct);
            D = [cA1' cA2' cB1' cB2' [0 0 0]'];
            
            %        D = INFLUENCECOEFF2(testDVE, DVECT, FP);
            
            velo(count,1:3) = [x_p, y_p, z_p];
            velo(count,4:6) = transpose(D*COEFF);
            
            count = count + 1;
            
        end
    end
end
hFig3 = figure(3);
clf(3)
patch(testDVE(:,1),testDVE(:,2),testDVE(:,3),'r','LineWidth',2)
alpha(0.5);
hold on
% scatter3(FP(1),FP(2),FP(3),100,'ok','filled');
quiver3(velo(:,1), velo(:,2), velo(:,3), velo(:,4), velo(:,5), velo(:,6));
hold off
grid on
axis equal
% axis tight
box on
xlabel('X-dir','FontSize',15);
ylabel('Y-dir','FontSize',15);
zlabel('Z-dir','FontSize',15);

%set proper zoom style
ax = gca;
z = zoom;
setAxes3DPanAndZoomStyle(z,ax,'camera')