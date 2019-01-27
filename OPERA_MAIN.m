clear
% clc
addpath('K Functions')

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
% strFILE = 'inputs/ellipse.dat';
% strFILE = 'inputs/box_wing.dat';
strFILE = 'inputs/goland_wing.dat';

[matPOINTS, strATYPE, vecSYM, flagRELAX, valMAXTIME, valDELTIME, valALPHA, valBETA, matTEPOINTS, matLEPOINTS, vecULS] = fcnOPREAD(strFILE);
[TR, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matELOC, matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG, matCONTROL, vecDVEAREA]...
    = fcnTRIANG(matPOINTS, 'SURFACE', []);
[vecLE, vecLEDVE, vecTE, vecTEDVE, matSPANDIR] = fcnLETEGEN(strATYPE, valNELE, matVLST, matELST, matDVECT, matEATT, matLEPOINTS, matTEPOINTS);

valALPHA = 30
valMAXTIME = 100
valDELTIME = 1
valDENSITY = 1.225;
flagRELAX =  0;

matUINF = repmat(fcnUINFWING(valALPHA, 0), valNELE, 1);

% [hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], [], matROTANG, [3 1 4 4], 'opengl');
% view([33, 28])
%% D-Matrix Creation
matKINCON_P = matCONTROL;
matKINCON_DVE = [1:valNELE]';

% Points where flow tangency is enforcedfe
matD = fcnDWING9(strATYPE, matEATT, matPLEX, valNELE, matELOC, matELST, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecLEDVE, vecTEDVE, vecSYM, matROTANG, matSPANDIR, matKINCON_P, matKINCON_DVE);
valDLEN = length(matD);

%% Preparing to timestep
% Initializing
matWELST = []; matWDVE = []; valWNELE = [];
matWEATT = []; matWEIDX = []; matWELOC = []; matWPLEX = []; matWDVECT = [];
matWALIGN = []; matWVATT = []; matWVNORM = []; matWCENTER = [];
matWAKEGEOM = []; matWCOEFF = []; matWVLST = []; matWROTANG = [];
vecWDVECIRC = []; CL = nan(valMAXTIME,1); CDi = nan(valMAXTIME,1);
e = nan(valMAXTIME,1); gust_vel_old = matCENTER.*0;

valWSIZE = length(vecTE);

% Building wing resultant
vecR = fcnRWING(valDLEN, 0, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matKINCON_P, matKINCON_DVE, matDVECT);

% Solving for wing coefficients
[matCOEFF] = fcnSOLVED(matD, vecR, valNELE);
%%
% Getting averaged mu values at vertices
% vmu = nan(size(matVLST,1),1);
% for i = 1:size(matVLST,1)
%     dves = matVATT(i, ~isnan(matVATT(i,:)))';
%     len = length(dves);
%     vloc = fcnGLOBSTAR(repmat(matVLST(i,:), len, 1) - matCENTER(dves,:), matROTANG(dves,:));
%     circ = sum([0.5.*vloc(:,2).^2 vloc(:,2) 0.5.*vloc(:,1).^2 vloc(:,1) ones(size(vloc(:,1)))].*matCOEFF(dves,:),2);
%     vmu(i) = mean(circ);
% end
% 
% % Getting averaged mu  values at edge midpoints
% emu = nan(size(matELST,1),1);
% for i = 1:size(matELST,1)
%     dves = matEATT(i, matEATT(i,:) > 0)';
%     len = length(dves);
%     pt = (matVLST(matELST(i,1),:) + matVLST(matELST(i,2),:))./2;
%     vloc = fcnGLOBSTAR(repmat(pt,len,1) - matCENTER(dves,:), matROTANG(dves,:));
%     circ = sum([0.5.*vloc(:,2).^2 vloc(:,2) 0.5.*vloc(:,1).^2 vloc(:,1) ones(size(vloc(:,1)))].*matCOEFF(dves,:),2);
%     emu(i) = mean(circ);
% end
% 
% % Set to zero leading edge and wing tips
% idx = ~all(matEATT,2);
% idx(vecTE) = 0;
% emu(idx) = 0;
% tmp = reshape(matELST(idx,:),[],1,1);
% vmu(tmp) = 0;
% 
% % Getting mu  values at DVE control points
% dmu = matCOEFF(:,5);


% [hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], [], matROTANG, [3 1 4 4], 'opengl');
% view([33, 28])
% hold on
% scatter3(matVLST(:,1), matVLST(:,2), -vmu, 'ok');
% scatter3(matCENTER(:,1), matCENTER(:,2), -dmu, '^m');
% pt = (matVLST(matELST(:,1),:) + matVLST(matELST(:,2),:))./2;
% scatter3(pt(:,1), pt(:,2), -emu, 'xr')
% hold off


% for i = 1:valNELE
%    vloc = fcnGLOBSTAR(matVLST(matDVE(i,:),:) - matCENTER(i,:), repmat(matROTANG(i,:), 3, 1));
%    vval = vmu(matDVE(i,:));
%    
%    eloc = (matVLST(matELST(matEIDX(i,:),1),:) + matVLST(matELST(matEIDX(i,:),2),:))./2;
%    eloc = fcnGLOBSTAR(eloc - matCENTER(i,:), repmat(matROTANG(i,:), 3, 1));
%    eval = emu(matEIDX(i,:));
%    
%    dloc = [0 0 0];
%    dval = matCOEFF(i,5);
%    
%    pts = [vloc; eloc; dloc]; 
%    XDATA = pts(:,1); YDATA = pts(:,2);
%    ZDATA = [vval; eval; dval]; 
%    
%    [xData, yData, zData] = prepareSurfaceData( XDATA, YDATA, ZDATA );
%    ft = fittype( 'poly22' );
%    [fitresult, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );
%    coeffvalues(fitresult);
% end

%%
matCOEFF_HSTRY(:,:,1) = matCOEFF;

valGUSTAMP = 0.0025;
valGUSTL = 7.3152;
flagGUSTMODE = 2;
valGUSTSTART = 12;

for valTIMESTEP = 1:valMAXTIME
    %% Timestep to solution
    %   Move wing
    %   Generate new wake elements
    %   Create and solve WD-Matrix for new elements
    %   Solve wing D-Matrix with wake-induced velocities
    %   Solve entire WD-Matrix (?)
    %   Relaxation procedure (Relax, create W-Matrix and W-Resultant, solve W-Matrix)
    %   Calculate surface normal forces
    %   Calculate DVE normal forces
    %   Calculate induced drag
    
    % Moving the wing
    [matVLST, matCENTER, matNEWWAKE, matCONTROL, matKINCON_P] = fcnMOVEWING(matUINF, valDELTIME, matVLST, matCENTER, matELST, vecTE, matCONTROL, matKINCON_P);
%     if valTIMESTEP >= valGUSTSTART
%     [matUINF, gust_vel_old] = fcnGUSTWING(matUINF, valGUSTAMP, valGUSTL, flagGUSTMODE, valDELTIME, 1, valGUSTSTART, matCENTER, gust_vel_old);
%     end
%     if valTIMESTEP == 10
%        strATYPE{3} = 'UNSTEADY';
%     end
    
    % Generating new wake elements
    if any(vecTE)
        [matWAKEGEOM, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC,...
            matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWCOEFF, matWROTANG, vecWLE, vecWTE, vecWLEDVE, vecWTEDVE, vecWDVECIRC] = fcnCREATEWAKE(valTIMESTEP, strATYPE, vecULS, matNEWWAKE, matWAKEGEOM, matCOEFF, valWSIZE, ...
            vecTE, vecTEDVE, matCENTER, matROTANG, matWCOEFF, matWPLEX, vecWDVECIRC);
        
        % Rebuild wing resultant
        vecR = fcnRWING(valDLEN, valTIMESTEP, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matKINCON_P, matKINCON_DVE, matDVECT);
        matCOEFF = fcnSOLVED(matD, vecR, valNELE);
        
        % Update wake coefficients
        if strcmpi(strATYPE{3},'STEADY')
            matWCOEFF = fcnDWAKE(strATYPE{3}, valTIMESTEP, strATYPE, vecULS, valWNELE, vecWLE, vecWLEDVE, vecWTE, vecWTEDVE, matWEATT, matWELST, matWROTANG, matWCENTER, matWVLST, vecTE, vecTEDVE, matCOEFF, matCENTER, matROTANG, matWCOEFF, matWPLEX, vecWDVECIRC);
        else
            [matWCOEFF(end - valWSIZE*2 + 1:end, :), vecWDVECIRC(end - valWSIZE*2 + 1:end, :)] = fcnDWAKENEW(valTIMESTEP, strATYPE, vecULS, valWNELE, vecWLE, vecWLEDVE, vecWTE, vecWTEDVE, matWEATT, matWELST, matWROTANG, matWCENTER, matWVLST, vecTE, vecTEDVE, matCOEFF, matCENTER, matROTANG, matWCOEFF, matWPLEX);
        end
        
        if flagRELAX == 1 && valTIMESTEP > 1
            [hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], matUINF, matROTANG, [3 1 4 4], 'opengl');
            [hFig1] = fcnPLOTWAKE(1, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER);
            setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera');
            view([-15 12])

            [matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG, matWAKEGEOM] = ...
                fcnRELAX2(valTIMESTEP, matWELOC, matWEATT, matWEIDX, matUINF, valDELTIME, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, ...
                matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST);
            
            %             [matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG, matWAKEGEOM] = ...
            %                 fcnRELAX(matUINF, valDELTIME, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, ...
            %                 matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST);
                        
            matWCOEFF = fcnDWAKE(valTIMESTEP, strATYPE, vecULS, valWNELE, vecWLE, vecWLEDVE, vecWTE, vecWTEDVE, matWEATT, matWELST, matWROTANG, matWCENTER, matWVLST, vecTE, vecTEDVE, matCOEFF, matCENTER, matROTANG, matWCOEFF, matWPLEX, vecWDVECIRC);
        end
        
        matCOEFF_HSTRY(:,:,valTIMESTEP + 1) = matCOEFF;
    end

    [CL(valTIMESTEP,1), CDi(valTIMESTEP,1), CY(valTIMESTEP,1), e(valTIMESTEP,1), vecDVELIFT, vecDVEDRAG] = fcnFORCES(valTIMESTEP, matVLST, matCENTER, matELST, matROTANG, ...
        matUINF, matCOEFF, vecTEDVE, valDENSITY, valNELE, matSPANDIR, vecTE, vecDVEAREA, matPLEX, matWCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG);
   
[hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], matUINF, matROTANG, [3 1 4 4], 'opengl');
[hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), matUINF, matROTANG, 'xr', 10);
% view([0 0])
% if valTIMESTEP > 0
%     [hFig1] = fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER);
%     [hFig1] = fcnPLOTCIRC(hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, matWPLEX, real(matWCOEFF), matUINF, matWROTANG, 'xb', 10);
% end
% view([-15 15])
    
%     fcnPLOTCOEFF(valTIMESTEP, matCOEFF_HSTRY);
%     hFig22 = figure(22);
%     clf(22);
%     set(hFig22,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
%     yyaxis right
%     plot(1:valTIMESTEP,CL(1:valTIMESTEP),'-^k');
%     ylabel('C_L','FontSize',15);
%     yyaxis left
%     plot(1:valTIMESTEP,CDi(1:valTIMESTEP),'--sk'); 
%     ylabel('C_D_i','FontSize',15);
%     xlabel('Timestep','FontSize',15);
%     grid minor
%     box on
%     legend('C_D_I','C_L','Location','SouthEast')
end
z
% Plot
[hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], matUINF, matROTANG, [3 1 4 4], 'opengl');
[hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), matUINF, matROTANG, 'xr', 10);
view([0 0])
if valTIMESTEP > 0
    [hFig1] = fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER);
    [hFig1] = fcnPLOTCIRC(hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, matWPLEX, real(matWCOEFF), matUINF, matWROTANG, 'xb', 10);
end
% setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera') ;

% s_ind = fcnSDVEVEL(matCENTER, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);
% w_ind = fcnSDVEVEL(matCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER);
% q_ind = s_ind + w_ind + matUINF;
% hold on
% quiver3(matCENTER(:,1), matCENTER(:,2), matCENTER(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3),'m')
% hold off

