clear
% clc
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
% strFILE = 'inputs/NACA 4412 2d.dat';
% strFILE = 'inputs/Circle_2d.dat';
strFILE = 'inputs/simple_wing.dat';
% strFILE = 'inputs/simple_wing2.dat';

[matPOINTS, strATYPE, vecSYM, flagRELAX, valMAXTIME, valDELTIME, valALPHA, valBETA, matTEPOINTS, matLEPOINTS] = fcnOPREAD(strFILE);

% matPOINTS = fcnSTLREAD('CAD Geom/master_airscrew.stl');
[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matELOC, matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG, matCONTROL, vecDVEAREA] = fcnTRIANG(matPOINTS);

flagRELAX = 0;
valMAXTIME = 0
valDELTIME = 0.1
valDENSITY = 1.225;

matUINF = repmat(fcnUINFWING(valALPHA, 0), valNELE, 1);

% [hFig1] = fcnPLOTBODY(1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], [], matROTANG, [3 1 4 4], 'opengl');
% view([33, 28])

%% D-Matrix Creation
vecTEDVE = [];
matSPANDIR = repmat([0 1 0],valNELE,1);
matSPANDIR = matSPANDIR - (dot(matSPANDIR, matDVECT(:,:,3),2)).*matDVECT(:,:,3);
vecLE = [];
vecTE = [];
vecLEDVE = [];
vecTEDVE = [];

if ~isempty(matLEPOINTS)
    [~, idxle(:,1)] = ismember(matLEPOINTS(:,:,1),matVLST,'rows');
    [~, idxle(:,2)] = ismember(matLEPOINTS(:,:,2),matVLST,'rows');
    [~, vecLE] = ismember(idxle, matELST,'rows');
    vecLEDVE = nonzeros(sort(matEATT(vecLE,:),2,'descend'));
end

matKINCON_P = matCONTROL;
matKINCON_DVE = [1:valNELE]';
if ~isempty(matTEPOINTS) && strcmpi(strATYPE{2},'PANEL')
    [~, idxte(:,1)] = ismember(matTEPOINTS(:,:,1),matVLST,'rows');
    [~, idxte(:,2)] = ismember(matTEPOINTS(:,:,2),matVLST,'rows');
    [~, vecTE] = ismember(idxte, matELST,'rows');
    vecTEDVE = nonzeros(sort(matEATT(vecTE,:),2,'descend'));
%     matKINCON_P = [matCONTROL; ((matVLST(matELST(vecTE,1),:) + matVLST(matELST(vecTE,2),:))./2)];
%     matKINCON_DVE = [[1:valNELE]'; matEATT(vecTE,1)];
elseif ~isempty(matTEPOINTS) && strcmpi(strATYPE{2},'THIN')
    [~, idxte(:,1)] = ismember(matTEPOINTS(:,:,1),matVLST,'rows');
    [~, idxte(:,2)] = ismember(matTEPOINTS(:,:,2),matVLST,'rows');
    [~, vecTE] = ismember(idxte, matELST,'rows');
    vecTEDVE = nonzeros(sort(matEATT(vecTE,:),2,'descend'));
%     matKINCON_P = [matCONTROL; ((matVLST(matELST(vecTE,1),:) + matVLST(matELST(vecTE,2),:))./2)];
%     matKINCON_DVE = [[1:valNELE]'; vecTEDVE];
else
    matKINCON_P = matCONTROL;
    matKINCON_DVE = [1:valNELE]';
end

% Points where flow tangency is enforced
matD = fcnDWING9(strATYPE, matEATT, matPLEX, valNELE, matELOC, matELST, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecLEDVE, vecSYM, matROTANG, matSPANDIR, matKINCON_P, matKINCON_DVE);

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
vecR = fcnRWING(valDLEN, 0, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matKINCON_P, matKINCON_DVE, matDVECT);

% Solving for wing coefficients
[matCOEFF] = fcnSOLVED(matD, vecR, valNELE);
matCOEFF_HSTRY(:,:,1) = matCOEFF;

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
    [matVLST, matCENTER, matNEWWAKE, matCONTROL] = fcnMOVEWING(matUINF, valDELTIME, matVLST, matCENTER, matELST, vecTE, matCONTROL);
    
    % Generating new wake elements
    if any(vecTE)
        [matWAKEGEOM, matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC,...
            matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWCOEFF, matWROTANG, vecWLE, vecWTE] = fcnCREATEWAKE(valTIMESTEP, matNEWWAKE, matWAKEGEOM, matCOEFF, valWSIZE, ...
            vecTE, vecTEDVE, matCENTER, matROTANG, matWCOEFF);
        
        % Rebuild wing resultant
        vecR = fcnRWING(valDLEN, valTIMESTEP, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matKINCON_P, matKINCON_DVE, matDVECT);
        matCOEFF = fcnSOLVED(matD, vecR, valNELE);
        
        matCOEFF_HSTRY(:,:,valTIMESTEP + 1) = matCOEFF;
    end
    
    if flagRELAX == 1
        [matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG] = ...
            fcnRELAX(matUINF, valDELTIME, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, ...
            matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST);
    end
    
    matDVENFORCE = fcnHDVENFORCE(strATYPE, matUINF, matCONTROL, matDVECT, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, vecDVEAREA);
    vecDVELIFT = sqrt(sum(matDVENFORCE.^2,2)).*cosd(valALPHA);
    vecDVEDRAG = sqrt(sum(matDVENFORCE.^2,2)).*sind(valALPHA);
    matDVELIFTDIR = cross(matUINF, matSPANDIR, 2);
    vecDVELIFT = dot( matDVELIFTDIR, matDVENFORCE, 2);
    vecDVEDRAG = dot(cross(matDVELIFTDIR, matSPANDIR, 2), matDVENFORCE, 2);
    CL = sum(vecDVELIFT)./(0.5.*sum(vecDVEAREA));
    CDi = sum(vecDVEDRAG)./(0.5.*sum(vecDVEAREA));
    fprintf('Timestep: %d\t\tCL = %0.5f\t\tCDi = %0.5f\n', valTIMESTEP, CL, CDi);
end

if valTIMESTEP > 0
    hFig21 = figure(21);
    clf(21);
    linestyles = {'--';'-.';'-';':'};
    markers = {'o';'x';'s';'^';'*';'d';'v';'>';'<';'p';'h'};
    colors = {'k';'b';'r';'m'};
    hold on
    for i = 1:5
        pFig1 = plot(1:valTIMESTEP,abs(reshape(sum(matCOEFF_HSTRY(:,i,2:end) - matCOEFF_HSTRY(:,i,1:end-1),1),[],1,1)),'-ok');
        pFig1.LineStyle = linestyles{1+mod(i,4),:};
        pFig1.Marker = markers{1+mod(i,11),:};
        pFig1.Color = colors{1+mod(i,4),:};
    end
    ylim([0, inf]);
    grid minor
    box on
    hold off
    xlabel('Timestep','FontSize',15)
    ylabel('Coefficient Delta','FontSize',15)
    legend('A_1','A_2','B_1','B_2','C_3','Location','NorthEast')
    set(gca, 'YScale', 'log')
end

%% Plot
[hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], matUINF, matROTANG, [3 1 4 4], 'opengl');
[hFig1, circ_all] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), matUINF, matROTANG, 'xr', 4);
if valTIMESTEP > 0
    [hFig1] = fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER);
        [hFig1,circ_all] = fcnPLOTCIRC(hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, matWPLEX, real(matWCOEFF), matUINF, matWROTANG, 'xb', 4);
end
% hold on
% quiver3(matCENTER(:,1), matCENTER(:,2), matCENTER(:,3), matDVENFORCE(:,1), matDVENFORCE(:,2), matDVENFORCE(:,3))
% hold off

% fpg = matCONTROL + 1e-6.*matDVECT(:,:,3);
% fpg = matCONTROL;

% q_inds = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);
% q_ind = q_inds + matUINF;

% fcolor = sqrt(sum(q_ind.^2,2));
% fcolor = 1 - (fcolor.^2);

% hold on
% p = patch('Faces',matDVE(:,:,1),'Vertices',matVLST,'FaceVertexCData',fcolor,'LineWidth',2);
% p.FaceColor = 'flat';
% quiver3(fpg(:,1),fpg(:,2),fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3), 'g')
% quiver3(fpg(:,1),fpg(:,2),fpg(:,3), q_inds(:,1), q_inds(:,2), q_inds(:,3), 'r')
% hold off
% colorbar;
% grid on
% box on
% axis equal
% hold off
% colorbar;
% view([33, 28])
% toc

granularity = .1;
x = -0.5:granularity:1.5;
% y = -3:granularity:3;
y = 0:granularity:3;
% y = 1.5;
z = -1:granularity:1;
[X,Y,Z] = meshgrid(x,y,z);
% y = x.*0 + 1.5;
y = [0.5 1.5 2.5];
fpg = [X(:) Y(:) Z(:)];
s_ind = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCONTROL);
q_ind = s_ind + repmat(matUINF(1,:), length(s_ind(:,1)),1);
Xq = reshape(q_ind(:,1), size(X));
Yq = reshape(q_ind(:,2), size(Y));
Zq = reshape(q_ind(:,3), size(Z));
% 
[Xs,Ys,Zs] = meshgrid(-0.5,y,z);
hold on
streamline(X,Y,Z,Xq,Yq,Zq,Xs,Ys,Zs);
% quiver3(fpg(:,1), fpg(:,2), fpg(:,3),q_ind(:,1),q_ind(:,2),q_ind(:,3));
hold off
% 
% % y = [0.5 1.5 2.5];
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
% x = 0.4916;
% y = 0.01667;
% z = 0.4:granularity:2;
% [X,Y,Z] = meshgrid(x,y,z);
% 
% fpg = unique([X(:), Y(:), Z(:)],'rows');
% s_ind = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);
% s_ind = s_ind + matUINF(1,:);
% 
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
% r = linspace(0.4, 2, 100);
% hold on
% scatter(sind(theta).*(1 + (0.5.^2./r.^2)), r, '^b')
% hold off
% 
% hFig21 = figure(22);
% clf(22);
% 
% scatter(s_ind(:,3), fpg(:,3), 'xk');
% grid minor
% box on
% axis tight
% hold on
% plot([min(s_ind(:,3)), max(s_ind(:,3))],[0.5 0.5],'-r')
% plot([min(s_ind(:,3)), max(s_ind(:,3))],[-0.5 -0.5],'-r')
% hold off
% xlabel('Normal Velocity','FontSize',15);
% ylabel('Z-Location','FontSize',15);
% 
% theta = 90;
% r = linspace(0.4, 2, 100);
% hold on
% scatter(cosd(theta).*(1 + (0.5.^2./r.^2)), r, '^b')
% hold off