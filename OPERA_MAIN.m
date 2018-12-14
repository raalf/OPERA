clear
clc
tic
%% Header
disp('====================================================================');
disp('                    /$$$$$$  /$$$$$$$  /$$$$$$$$ /$$$$$$$   /$$$$$$ ');
disp('+---------------+  /$$__  $$| $$__  $$| $$_____/| $$__  $$ /$$__  $$');
disp('| RYERSON       | | $$  \ $$| $$  \ $$| $$      | $$  \ $$| $$  \ $$');
disp('| APPLIED       | | $$  | $$| $$$$$$$/| $$$$$   | $$$$$$$/| $$$$$$$$');
disp('| AERODYNAMICS  | | $$  | $$| $$____/ | $$__/   | $$__  $$| $$__  $$');
disp('| LABORATORY OF | | $$  | $$| $$      | $$      | $$  \ $$| $$  | $$');
disp('| FLIGHT        | |  $$$$$$/| $$      | $$$$$$$$| $$  | $$| $$  | $$');
disp('+---------------+  \______/ |__/      |________/|__/  |__/|__/  |__/');
disp('====================================================================');
%% Preamble
% strFILE = 'inputs/simple_wing2d.dat';
% strFILE = 'inputs/2dve.dat';
strFILE = 'inputs/NACA 4412 2d.dat';
% strFILE = 'inputs/Circle_2d.dat';

[matPOINTS, strATYPE, vecSYM, flagRELAX, valMAXTIME, valDELTIME, valALPHA, valBETA, matTEPOINTS, matLEPOINTS] = fcnOPREAD(strFILE);

% matPOINTS = fcnSTLREAD('CAD Geom/master_airscrew.stl');
[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matELOC, matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG, matCONTROL] = fcnTRIANG(matPOINTS);

flagRELAX = 0;
valMAXTIME = 0;
valDENSITY = 1.225;

matUINF = repmat(fcnUINFWING(valALPHA, 0), valNELE, 1);

% [hFig1] = fcnPLOTBODY(1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], [], matROTANG, [3 1 4 4], 'opengl');
% view([33, 28])

%% D-Matrix Creation
vecTEDVE = [];
vecSPANDIR = repmat([0 1 0],valNELE,1);
vecSPANDIR = vecSPANDIR - (dot(vecSPANDIR, matDVECT(:,:,3),2)).*matDVECT(:,:,3);
vecLE = [1];
vecTE = [];
vecLEDVE = nonzeros(sort(matEATT(vecLE,:),2,'descend'));
if ~isempty(vecTE)
    vecTEDVE = nonzeros(sort(matEATT(vecTE,:),2,'descend')); % A vector of trailing edge HDVEs, which corresponds to vecTE edges
    vecSPANDIR = fcnGLOBSTAR(repmat([0 1 0], length(vecTEDVE)), matROTANG(vecTEDVE,:)); % Spanwise direction for each HDVE (may change with rotor stuff)
end

matD = fcnDWING9(strATYPE, matEATT, matPLEX, valNELE, matELOC, matELST, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecLEDVE, vecSYM, matROTANG, vecSPANDIR, vecTEDVE, matCONTROL);

valDLEN = length(matD);

%% Preparing to timestep
matWADJE = [];
matWELST = [];
matWDVE = [];
valWNELE = [];
matWEATT = [];
matWEIDX = [];
matWELOC = [];
matWPLEX = [];
matWDVECT = [];
matWALIGN = [];
matWVATT = [];
matWVNORM = [];
matWCENTER = [];
matWAKEGEOM = [];
valWSIZE = [];
matWCOEFF = [];
matWVLST = [];
matWROTANG = [];

% Building wing resultant
vecR = fcnRWING(valDLEN, 0, matCENTER, matDVECT, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER);

% Solving for wing coefficients
[matCOEFF] = fcnSOLVED(matD, vecR, valNELE);

for valTIMESTEP = 1:valMAXTIME
    %% Timestep to solution
    %   Move wing
    %   Generate new wake elements
    %   Create and solve WD-Matrix for new elements
    %   Solve wing D-Matrix with wake-induced velocities
    %   Solve entire WD-Matrix
    %   Relaxation procedure (Relax, create W-Matrix and W-Resultant, solve W-Matrix)
    %   Calculate surface normal forces
    %   Calculate DVE normal forces
    %   Calculate induced drag
    
    valWSIZE = length(vecTE);
    
    % Moving the wing
    [matVLST, matCENTER, matNEWWAKE] = fcnMOVEWING(matUINF, valDELTIME, matVLST, matCENTER, matELST, vecTE);
    
    % Generating new wake elements
    if any(vecTE)
        [matWAKEGEOM, matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC,...
            matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWCOEFF, matWROTANG] = fcnCREATEWAKE(valTIMESTEP, matNEWWAKE, matWAKEGEOM, matCOEFF, valWSIZE, ...
            vecTE, vecTEDVE, matPLEX, matELOC, matELST, matDVE, matWCOEFF);
        
        % Rebuild wing resultant
        vecR = fcnRWING(valDLEN, valTIMESTEP, matCENTER, matDVECT, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER);
        matCOEFF = fcnSOLVED(matD, vecR, valNELE);
       
    end
         
end

%% Plot
[hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], matUINF, matROTANG, [3 1 4 4], 'opengl');
[hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), matUINF, matROTANG, 'r', 1);

fpg = matCENTER + 0.0001.*matDVECT(:,:,3);
% fpg = matCENTER;

q_inds = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);
q_ind = q_inds + matUINF;
% q_inds = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCONTROL);
% q_ind = -fcnSTARGLOB([zeros(size(matCOEFF(:,2),1),1) matCOEFF(:,2)./8 zeros(size(matCOEFF(:,2),1),1)], matROTANG);

fcolor = sqrt(sum(q_ind.^2,2));
fcolor = 1 - (fcolor.^2);

% fcolor = [matCOEFF(:,2)./2 zeros(size(matCOEFF(:,2),1),2)] + matUINF;
% fcolor = sqrt(sum(fcolor.^2,2));
% velocities = [matCOEFF(:,2), matCOEFF(:,4), matCOEFF(:,1).*0];
% fcolor = 1 - (sqrt(sum(velocities.^2,2)) + 1).^2;

hold on
p = patch('Faces',matDVE(:,:,1),'Vertices',matVLST,'FaceVertexCData',fcolor,'LineWidth',2);
p.FaceColor = 'flat';
quiver3(fpg(:,1),fpg(:,2),fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3), 'g')
quiver3(fpg(:,1),fpg(:,2),fpg(:,3), q_inds(:,1), q_inds(:,2), q_inds(:,3), 'r')
hold off
colorbar;
grid on
box on
axis equal
hold off
colorbar;
if valTIMESTEP > 0
    [hFig1] = fcnPLOTWAKE(1, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER);
    [hFig1] = fcnPLOTCIRC(hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, matWPLEX, real(matWCOEFF), vecUINF, matWROTANG, 'b', 20);
end
view([33, 28])
toc

% granularity = .05;
% x = -0.5:granularity:1.5;
% % y = -3:granularity:3;
% y = -1.5:granularity:1.5;
% % y = 0;
% z = -1:granularity:1;
% [X,Y,Z] = meshgrid(x,y,z);
% y = x.*0 + 0.025;
% 
% s_ind = fcnSDVEVEL([X(:) Y(:) Z(:)], valNELE, matCOEFF, matPLEX, matROTANG, matCONTROL);
% q_ind = s_ind + repmat(matUINF(1,:), length(s_ind(:,1)),1);
% Xq = reshape(q_ind(:,1), size(X));
% Yq = reshape(q_ind(:,2), size(Y));
% Zq = reshape(q_ind(:,3), size(Z));
% 
% [Xs,Ys,Zs] = meshgrid(-0.5,y,z);
% hold on
% streamline(X,Y,Z,Xq,Yq,Zq,Xs,Ys,Zs);
% % quiver3(X(:),Y(:),Z(:),q_ind(:,1),q_ind(:,2),q_ind(:,3));
% hold off

% hFig20 = figure(20);
% clf(20);
% scatter(matCENTER(:,1), fcolor,'xr')
% hold on
% theta = linspace(0, pi, 100);
% plot(((1 - cos(theta)))./2, 1 - 4.*sin(theta).^2,'--ok')
% hold off
% set(gca,'Ydir','reverse')
% grid minor
% box on
% axis tight

% granularity = .01;
% x = 0.5;
% y = 0.025;
% z = -1:granularity:1;
% [X,Y,Z] = meshgrid(x,y,z);
% 
% fpg = unique([X(:), Y(:), Z(:)],'rows');
% s_ind = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);
% s_ind = s_ind + matUINF(1,:);
% hFig21 = figure(21);
% clf(21);
% 
% scatter(s_ind(:,1), fpg(:,3), 'xk');
% grid minor
% box on
% axis tight
% hold on
% plot([min(s_ind(:,1)), max(s_ind(:,1))],[0.5 0.5],'-r')
% plot([min(s_ind(:,1)), max(s_ind(:,1))],[-0.5 -0.5],'-r')
% hold off
% xlabel('Tangential Velocity','FontSize',15);
% ylabel('Z-Location','FontSize',15);
% 
% theta = 90;
% r = linspace(0.2, 1, 100);
% hold on
% scatter(sind(theta).*(1 + (0.5.^2./r.^2)), r, '^b')
% hold off
