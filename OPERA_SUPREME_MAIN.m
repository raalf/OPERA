% clear

%% Header
disp('=========================================================================');
disp('                         /$$$$$$  /$$$$$$$  /$$$$$$$$ /$$$$$$$   /$$$$$$ ');
disp('  +---------------+     /$$__  $$| $$__  $$| $$_____/| $$__  $$ /$$__  $$');
disp('  | RYERSON       |    | $$  \ $$| $$  \ $$| $$      | $$  \ $$| $$  \ $$');
disp('  | APPLIED       |    | $$  | $$| $$$$$$$/| $$$$$   | $$$$$$$/| $$$$$$$$');
disp('  | AERODYNAMICS  |    | $$  | $$| $$____/ | $$__/   | $$__  $$| $$__  $$');
disp('  | LABORATORY OF |    | $$  | $$| $$      | $$      | $$  \ $$| $$  | $$');
disp('  | FLIGHT        |    |  $$$$$$/| $$      | $$$$$$$$| $$  | $$| $$  | $$');
disp('  +---------------+     \______/ |__/      |________/|__/  |__/|__/  |__/');
disp('  /$$$$$$  /$$   /$$ /$$$$$$$  /$$$$$$$  /$$$$$$$$ /$$      /$$ /$$$$$$$$ ');
disp(' /$$__  $$| $$  | $$| $$__  $$| $$__  $$| $$_____/| $$$    /$$$| $$_____/ ');
disp('| $$  \__/| $$  | $$| $$  \ $$| $$  \ $$| $$      | $$$$  /$$$$| $$       ');
disp('|  $$$$$$ | $$  | $$| $$$$$$$/| $$$$$$$/| $$$$$   | $$ $$/$$ $$| $$$$$    ');
disp(' \____  $$| $$  | $$| $$____/ | $$__  $$| $$__/   | $$  $$$| $$| $$__/    ');
disp(' /$$  \ $$| $$  | $$| $$      | $$  \ $$| $$      | $$\  $ | $$| $$      ');
disp('|  $$$$$$/|  $$$$$$/| $$      | $$  | $$| $$$$$$$$| $$ \/  | $$| $$$$$$$$ ');
disp(' \______/  \______/ |__/      |__/  |__/|________/|__/     |__/|________/ ');
disp('=========================================================================');

%% Reading in geometry
filename = 'inputs/TMotor_Coarse2_offset.vap'
% filename = 'inputs/QuadRotor.vap'
% filename = 'inputs/kussner.vap'

[flgRELAX, flgSTEADY, valMAXTIME, valDELTIME, valDENSITY, valUINF, valALPHA, valBETA, ...
    valROLL, valFPA, valTRACK, valAREA, valSPAN, matPOINTS, matTEPOINTS, matLEPOINTS, ...
    vecDVESURFACE, vecDVEFLIP, valROTORS, valWINGS, vecROTORRPM, vecROTORDIAM, ...
    vecROTORBLADES, matROTORHUB, matROTORAXIS, vecDVEWING, vecDVEROTOR, vecDVESDFLIP, vecVEHORIG] = fcnXMLREAD(filename);

[TR, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matPLEX, matDVECT, matVATT, ~, ...
    matCENTER, matROTANG, ~, vecDVEAREA, matSPANDIR, vecDVECHORD] = fcnTRIANG(matPOINTS, vecDVEFLIP, vecDVESDFLIP);

flgGIF = false;

valALPHA = SEQ_ALPHAS(jjj)
valUINF = SEQ_VS(jjjj)

% Rolling
% valROTSTART = 0 % Starting timestep
% vecROTRATE = [50 0 0] % [roll, pitch, yaw] rad/s
valROTSTART = 120; % Starting timestep
vecROTRATE = [SEQ_RATES(jjjjj) 0 0] % [roll, pitch, yaw] rad/s

% valGUSTAMP = 0.053.*valUINF;
% valGUSTL = 100;
% valGUSTMODE = 3;
% valGUSTSTART = 1;

% disp('HUB TRANSLATION APPROXIMATION EVALUATION in fcnUINF')

%% Preliminary Geometry Stuff
% Duplicating rotor blades if need be
if valROTORS > 0
    [valNELE, matVLST, matELST, matDVE, matCENTER, matPLEX, vecDVEAREA, ...
        matEATT, matEIDX, vecDVESURFACE, vecDVEFLIP, matDVECT, ...
        matROTANG, matVATT, matTEPOINTS, matLEPOINTS, matSPANDIR, vecDVEROTOR, vecDVEWING, vecDVESDFLIP, matROTORTRANS, vecDVECHORD] ...
        = fcnREADYROTOR(valROTORS, valNELE, matVLST, matELST, matDVE, ...
        vecDVEAREA, matEATT, matEIDX, vecDVESURFACE, vecDVEFLIP, matVATT, ...
        matTEPOINTS, matLEPOINTS, matSPANDIR, vecROTORBLADES, matROTORHUB, matROTORAXIS, ...
        vecDVEWING, vecDVEROTOR, vecDVESDFLIP, vecDVECHORD);
    vecJ = abs(valUINF)./((vecROTORRPM.*(pi/30)).*(vecROTORDIAM/2));
else
    matROTORTRANS = [];
end
% hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], [], matROTANG, [], 'opengl');

% Generating leading and trailing edge information
[vecLE, vecLEDVE, vecTE, vecTEDVE, matELST, vecCHORD, vecSPAN] = fcnLETEGEN(matVLST, matELST, matEATT, matLEPOINTS, matTEPOINTS);

% Kinematic condition locations
[matKINCON_P, vecKINCON_DVE] = fcnKINCON(valNELE, matDVE, matCENTER, matVLST, vecDVEAREA);

% DVE velocities
[matUINF, matVUINF, matUINF_KK] = fcnUINF(valUINF, valALPHA, valBETA, valNELE, valROTORS, matVLST, matKINCON_P, ...
    vecROTORRPM, vecDVEROTOR, matCENTER, matDVE, matROTORHUB, matROTORAXIS, vecKINCON_DVE, valROTSTART, vecROTRATE, 0, vecVEHORIG);

% Force directions
matDRAG_DIR = matUINF./sqrt(sum(matUINF(1,:).^2,2));
matLIFT_DIR = cross(matDRAG_DIR, matSPANDIR, 2);
matLIFT_DIR = matLIFT_DIR./sqrt(sum(matLIFT_DIR.^2,2));
matSIDE_DIR = cross(matLIFT_DIR, matDRAG_DIR, 2);
matSIDE_DIR = matSIDE_DIR./sqrt(sum(matSIDE_DIR.^2,2));

%% Creating GIF file
if flgGIF == true
    if exist('GIF','dir') ~= 7
        mkdir('GIF');
    end
    clear dir
    valGIFNUM = size(dir('GIF'),1) - 1;
end

%% D-Matrix Creation
matD = fcnDWING9(matEATT, matPLEX, valNELE, matELST, matVLST, matCENTER, matDVECT, vecTE, vecLE, vecLEDVE, vecTEDVE, matROTANG, matKINCON_P, vecKINCON_DVE, vecDVESDFLIP);
valDLEN = length(matD);

%% Preparing to timestep
valZTOL = 1e-5; % Used to offset field points from surface sheets when necessary

% Initializing
matWELST = []; matWDVE = []; valWNELE = []; matWEATT = []; matWEIDX = []; matWPLEX = []; matWDVECT = [];
matWCENTER = []; matWAKEGEOM = []; matWCOEFF = []; matWVLST = []; matWROTANG = []; vecWVMU = []; vecWEMU = [];
matWVGRID = []; matWEGRID = []; vecWDVEFLIP = logical([]); vecWLEDVE = []; vecWTEDVE = [];
CL = nan(valMAXTIME,1); CDi = nan(valMAXTIME,1); CT = nan(valMAXTIME,valROTORS); CMx = nan(valMAXTIME,valROTORS);
CMy = nan(valMAXTIME,valROTORS); CMz = nan(valMAXTIME,valROTORS); e = nan(valMAXTIME,1);
valWSIZE = length(vecTE); vecWDVESURFACE = []; matINTCIRC = nan(valNELE,valMAXTIME); vecTSITER = nan(valMAXTIME,2);
vecRSQUARED = nan(valMAXTIME,2); matWDVEGRID = []; vecWDVESDFLIP = logical([]); gust_vel_old = zeros(size(vecKINCON_DVE));

% Building wing resultant
vecR = fcnRWING(flgSTEADY, valZTOL, valDLEN, 0, matUINF_KK, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, matKINCON_P, vecKINCON_DVE, matDVECT, [], []);
boolKINCON = false(size(matD,1),1);
boolKINCON(end-(valNELE*4-1):end) = true;

% Solving for wing coefficients
matCOEFF = fcnSOLVED(matD, vecR, valNELE);
[vecVMU, vecEMU] = fcnVEMU(matVLST, matVATT, matCENTER, matROTANG, matCOEFF, matELST, matEATT, vecTE);

%% Timestep to solution
for valTIMESTEP = 1:valMAXTIME
    tic
    [matVLST, matCENTER, matNEWWAKE, matUINF, matVUINF, matUINF_KK, matPLEX, matDVECT, matROTANG, matSPANDIR, matKINCON_P, matROTORHUB, matROTORAXIS, matROTORTRANS, vecVEHORIG] = ...
        fcnMOVESURFACES(valNELE, valUINF, valROTORS, valALPHA, valBETA, vecROTORRPM, valDELTIME, matROTORHUB, matVLST, matELST, vecTE, matDVE, ...
        vecDVEFLIP, matSPANDIR, matKINCON_P, vecDVEROTOR, vecKINCON_DVE, matROTORAXIS, valROTSTART, vecROTRATE, valTIMESTEP, matROTORTRANS, vecVEHORIG);
%     
%         if ~isempty(valGUSTSTART) && valGUSTSTART >= 0
%             [matUINF_KK, gust_vel_old] = fcnGUSTWING(matUINF_KK, valGUSTAMP, valGUSTL, valGUSTMODE, valDELTIME, valUINF, valGUSTSTART, matKINCON_P, gust_vel_old, 0);
%         end
    
    % Generating new wake elements
    if any(vecTE)
        [matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWPLEX, ...
            matWDVECT, matWCENTER, matWCOEFF, matWROTANG, vecWLE, vecWTE, vecWLEDVE, vecWTEDVE, ...
            matWVGRID, vecWOTE, vecWOTEDVE, matWEGRID, matWE2GRID, vecWVMU, vecWEMU, vecWDVEFLIP, ...
            vecWDVESURFACE, matWDVEGRID, vecWDVESDFLIP] = fcnCREATEWAKE2(valTIMESTEP, flgSTEADY, matNEWWAKE, matCOEFF, valWSIZE, ...
            vecTEDVE, matCENTER, matROTANG, matWCOEFF, matWPLEX, matWVGRID, matWVLST, matWELST, matWEGRID, ...
            vecWVMU, vecWEMU, vecWDVEFLIP, matWCENTER, matWROTANG, matWDVECT, ...
            matWDVE, matWEIDX, vecWLEDVE, matWEATT, vecDVESURFACE, vecWDVESURFACE, matWDVEGRID, vecWDVESDFLIP, vecDVESDFLIP);
        
        %         matD = fcnDWING9(matEATT, matPLEX, valNELE, matELST, matVLST, matCENTER, matDVECT, vecTE, vecLE, vecLEDVE, vecTEDVE, matROTANG, matKINCON_P, vecKINCON_DVE, vecDVESDFLIP);
        
        [vecR, matCOEFF, vecVMU, vecEMU, matWCOEFF, vecWVMU, vecWEMU, vecTSITER(valTIMESTEP,1), vecRSQUARED] = ...
            fcnTSITER(matCOEFF, valDLEN, valTIMESTEP, matUINF_KK, valWNELE, ...
            matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matKINCON_P, vecKINCON_DVE, ...
            matDVECT, matD, vecR, valNELE, matVLST, matVATT, matCENTER, matROTANG, ...
            matDVE, matELST, matEIDX, flgSTEADY, vecWLE, matWVGRID, matWEGRID, matWE2GRID, ...
            vecWVMU, vecWEMU, matWELST, matWVLST, vecTEDVE, vecWOTE, matWDVE, matWEIDX, ...
            matEATT, boolKINCON, vecTE, vecVMU, vecEMU, valZTOL, vecRSQUARED, vecWDVESDFLIP);
        
        % Relaxing Wake
        if flgRELAX == true && valTIMESTEP > 2
            [matWELST, matWVLST, matWPLEX, matWDVECT, matWCENTER, matWROTANG] = ...
                fcnRELAX8(valTIMESTEP, valDELTIME, valNELE, matCOEFF, matPLEX, valWNELE, matWCOEFF, matWDVE, matWVLST, matWPLEX, ...
                matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, matWELST, matWVGRID, vecWDVEFLIP, matWDVEGRID, vecWDVESDFLIP, vecDVESDFLIP);
            
            matWCOEFF = fcnADJCOEFF(vecWVMU, vecWEMU, matWVLST, matWCENTER, matWROTANG, matWDVE, matWCOEFF, matWELST, matWEIDX, valWNELE);
            
            [vecR, matCOEFF, vecVMU, vecEMU, matWCOEFF, vecWVMU, vecWEMU, vecTSITER(valTIMESTEP,2), vecRSQUARED] = ...
                fcnTSITER(matCOEFF, valDLEN, valTIMESTEP, matUINF_KK, valWNELE, ...
                matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matKINCON_P, vecKINCON_DVE, ...
                matDVECT, matD, vecR, valNELE, matVLST, matVATT, matCENTER, matROTANG, ...
                matDVE, matELST, matEIDX, flgSTEADY, vecWLE, matWVGRID, matWEGRID, matWE2GRID, ...
                vecWVMU, vecWEMU, matWELST, matWVLST, vecTEDVE, vecWOTE, matWDVE, matWEIDX, ...
                matEATT, boolKINCON, vecTE, vecVMU, vecEMU, valZTOL, vecRSQUARED, vecWDVESDFLIP);
        end
        
        % Updating coefficient convergence history
        matCOEFF_HSTRY(:,:,valTIMESTEP) = matCOEFF;
        
        if flgGIF == true
            hFig1 = fcnGIF(valTIMESTEP, valNELE, matDVE, matVLST, matCENTER, matELST, matDVECT, matPLEX, matCOEFF, matUINF, matROTANG, ...
                valWNELE, matWDVE, matWVLST, matWCENTER, matWELST, matWDVECT, vecWDVESURFACE, valWSIZE, matWVGRID, valGIFNUM);
        end      
                
        %% Calculating Forces
        [matF_FS(:,:,valTIMESTEP), matF_IF(:,:,valTIMESTEP), matF_ID(:,:,valTIMESTEP), matF_NC(:,:,valTIMESTEP), matINTCIRC(:,valTIMESTEP), matKINCON_V_ALL(:,:,valTIMESTEP)] = ...
            fcnDVEFORCES2(matVLST, matELST, matROTANG, ...
            matCOEFF, vecTEDVE, valDENSITY, valNELE, matSPANDIR, vecTE, matPLEX, matWCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG, ...
            matWVLST, vecWLEDVE, matWELST, valWSIZE, matWDVE, vecWDVEFLIP, matWEIDX, vecWEMU, ...
            vecWVMU, matCENTER, matKINCON_P, vecKINCON_DVE, matUINF_KK, matDVECT, valZTOL, vecWDVESDFLIP, vecDVESDFLIP, valDELTIME, matCOEFF_HSTRY, valTIMESTEP, vecDVECHORD);
        
        [CL(valTIMESTEP), CDi(valTIMESTEP), e(valTIMESTEP), ~, ~, ~, ~, ~, CT(valTIMESTEP,:), CMx(valTIMESTEP,:), CMy(valTIMESTEP,:), CMz(valTIMESTEP,:)] = ...
            fcnFORCES(true, valTIMESTEP, valWINGS, valROTORS, ...
            matF_FS, matF_IF, matF_ID, [], matLIFT_DIR, matSIDE_DIR, matDRAG_DIR, valDENSITY, valAREA, valSPAN, valUINF, vecROTORRPM, ...
            vecROTORDIAM, matROTORAXIS, matROTORHUB, matCENTER, vecTSITER, vecRSQUARED, vecDVEROTOR, vecDVEWING, matROTORTRANS);
                
        matSPANDIR_ALL(:,:,valTIMESTEP) = matSPANDIR;
        matUINF_ALL(:,:,valTIMESTEP) = matUINF;
        matUINF_KK_ALL(:,:,valTIMESTEP) = matUINF_KK;
        matKINCON_P_ALL(:,:,valTIMESTEP) = matKINCON_P;
        
        matROTANG_ALL(:,:,valTIMESTEP) = matROTANG;
        matCENTER_ALL(:,:,valTIMESTEP) = matCENTER;
        matROTORAXIS_ALL(:,:,valTIMESTEP) = matROTORAXIS;
        matROTORHUB_ALL(:,:,valTIMESTEP) = matROTORHUB;
        matROTORTRANS_ALL(:,:,:,valTIMESTEP) = matROTORTRANS;
        
    end
    
    vecTSTIME(valTIMESTEP) = toc;
%     save(['timesteps/timestep_', num2str(valTIMESTEP), '.mat']);
end

%% Apparent mass
skip = 1;
matF_NCP = fcnDGAMMADT(skip, valTIMESTEP, valNELE, valDENSITY, valDELTIME, matINTCIRC, matDVECT);

for ii = 1:skip:valMAXTIME
    [CL_UP(ii,1), ~, ~, ~, ~, ~, ~, ~, CT_UP(ii,:), CMx_UP(ii,:), CMy_UP(ii,:), CMz_UP(ii,:)] = fcnFORCES(false, ii, valWINGS, valROTORS, ...
        matF_FS, matF_IF, matF_ID, matF_NCP, matLIFT_DIR, matSIDE_DIR, matDRAG_DIR, valDENSITY, valAREA, valSPAN, valUINF, vecROTORRPM, ...
        vecROTORDIAM, matROTORAXIS_ALL(:,:,ii), matROTORHUB_ALL(:,:,ii), matCENTER_ALL(:,:,ii), ...
        vecTSITER, vecRSQUARED, vecDVEROTOR, vecDVEWING, matROTORTRANS_ALL(:,:,:,ii));
end

%% Save environment
save(['run_',num2str(now),'.mat']);

%% Post
% hFig2 = figure(2);
% clf(2);
% plot(CMx)
% hold on; plot((1:skip:valMAXTIME), CMx_UP(1:skip:valMAXTIME), '--b'); hold off
% 
% hFig3 = figure(3);
% clf(3);
% plot(CMy)
% hold on; plot((1:skip:valMAXTIME), CMy_UP(1:skip:valMAXTIME), '--b'); hold off
% 
% hFig4 = figure(4);
% clf(4);
% plot(CT,'-k')
% hold on; plot(CT_UP, '--b'); hold off
