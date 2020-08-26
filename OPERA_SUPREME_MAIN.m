clear

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
filename = 'inputs/TMotor_Coarse2.vap'
% filename = 'inputs/WIPP_FINAL_Coarse.vap'

[flgRELAX, flgSTEADY, valMAXTIME, valDELTIME, valDENSITY, valUINF, valALPHA, valBETA, ...
    valROLL, valFPA, valTRACK, valAREA, valSPAN, matPOINTS, matTEPOINTS, matLEPOINTS, ...
    vecDVESURFACE, vecDVEFLIP, valROTORS, valWINGS, vecROTORRPM, vecROTORDIAM, ...
    vecROTORBLADES, matROTORHUB, matROTORAXIS, vecDVEWING, vecDVEROTOR, vecDVESDFLIP] = fcnXMLREAD(filename);

[TR, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matPLEX, matDVECT, matVATT, ~, ...
    matCENTER, matROTANG, ~, vecDVEAREA, matSPANDIR, vecDVECHORD] = fcnTRIANG(matPOINTS, vecDVEFLIP, vecDVESDFLIP);

flgGIF = false;

% valUINF = Vs(jjj)
% vecROTORRPM = sign(vecROTORRPM).*repmat(RPMs(jjj), length(vecROTORRPM), 1)
% valALPHA = Alphas(jjjj)

%% Preliminary Geometry Stuff
% Duplicating rotor blades if need be
if valROTORS > 0
    [valNELE, matVLST, matELST, matDVE, matCENTER, matPLEX, vecDVEAREA, ...
        matEATT, matEIDX, vecDVESURFACE, vecDVEFLIP, matDVECT, ...
        matROTANG, matVATT, matTEPOINTS, matLEPOINTS, matSPANDIR, vecDVEROTOR, vecDVEWING, vecDVESDFLIP] ...
        = fcnREADYROTOR(valROTORS, valNELE, matVLST, matELST, matDVE, ...
        vecDVEAREA, matEATT, matEIDX, vecDVESURFACE, vecDVEFLIP, matVATT, ...
        matTEPOINTS, matLEPOINTS, matSPANDIR, vecROTORBLADES, matROTORHUB, matROTORAXIS, ...
        vecDVEWING, vecDVEROTOR, vecDVESDFLIP);
end
% hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], [], matROTANG, [], 'opengl');

% Generating leading and trailing edge information
[vecLE, vecLEDVE, vecTE, vecTEDVE, matELST, vecCHORD, vecSPAN] = fcnLETEGEN(matVLST, matELST, matEATT, matLEPOINTS, matTEPOINTS);

% Kinematic condition locations
[matKINCON_P, vecKINCON_DVE] = fcnKINCON(valNELE, matDVE, matCENTER, matVLST, vecDVEAREA);

% DVE velocities
[matUINF, matVUINF, matUINF_KK] = fcnUINF(valUINF, valALPHA, valBETA, valNELE, valROTORS, matVLST, matKINCON_P, ...
    vecROTORRPM, vecDVEROTOR, matCENTER, matDVE, matROTORHUB, matROTORAXIS, vecKINCON_DVE);

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
CL = nan(valMAXTIME,1); CDi = nan(valMAXTIME,1); CT = nan(valMAXTIME,valROTORS); e = nan(valMAXTIME,1);
valWSIZE = length(vecTE); vecWDVESURFACE = []; matINTCIRC = nan(valNELE,valMAXTIME); vecTSITER = nan(valMAXTIME,2);
vecRSQUARED = nan(valMAXTIME,2); matWDVEGRID = []; vecWDVESDFLIP = logical([]);

% Building wing resultant
vecR = fcnRWING(flgSTEADY, valZTOL, valDLEN, 0, matUINF_KK, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER, matKINCON_P, vecKINCON_DVE, matDVECT, [], []);
boolKINCON = false(size(matD,1),1);
boolKINCON(end-(valNELE*4-1):end) = true;

% Solving for wing coefficients
matCOEFF = fcnSOLVED(matD, vecR, valNELE);
[vecVMU, vecEMU] = fcnVEMU(matVLST, matVATT, matCENTER, matROTANG, matCOEFF, matELST, matEATT, vecTE);
matCOEFF_HSTRY(:,:,1) = matCOEFF;


%% Timestep to solution
for valTIMESTEP = 1:valMAXTIME
    tic
    
    [matVLST, matCENTER, matNEWWAKE, matUINF, matVUINF, matUINF_KK, matPLEX, matDVECT, matROTANG, matSPANDIR, matKINCON_P, matROTORHUB] = ...
        fcnMOVESURFACES(valNELE, valUINF, valROTORS, valALPHA, valBETA, vecROTORRPM, valDELTIME, matROTORHUB, matVLST, matELST, vecTE, matDVE, ...
        vecDVEFLIP, matSPANDIR, matKINCON_P, vecDVEROTOR, vecKINCON_DVE, matROTORAXIS);
    
    % Generating new wake elements
    if any(vecTE)
        [matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWPLEX, ...
            matWDVECT, matWCENTER, matWCOEFF, matWROTANG, vecWLE, vecWTE, vecWLEDVE, vecWTEDVE, ...
            matWVGRID, vecWOTE, vecWOTEDVE, matWEGRID, matWE2GRID, vecWVMU, vecWEMU, vecWDVEFLIP, ...
            vecWDVESURFACE, matWDVEGRID, vecWDVESDFLIP] = fcnCREATEWAKE2(valTIMESTEP, flgSTEADY, matNEWWAKE, matCOEFF, valWSIZE, ...
            vecTEDVE, matCENTER, matROTANG, matWCOEFF, matWPLEX, matWVGRID, matWVLST, matWELST, matWEGRID, ...
            vecWVMU, vecWEMU, vecWDVEFLIP, matWCENTER, matWROTANG, matWDVECT, ...
            matWDVE, matWEIDX, vecWLEDVE, matWEATT, vecDVESURFACE, vecWDVESURFACE, matWDVEGRID, vecWDVESDFLIP, vecDVESDFLIP);
        
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
                matEATT, boolKINCON, vecTE, vecVMU, vecEMU, valZTOL, vecRSQUARED);
        end
        
        % Updating coefficient convergence history
        matCOEFF_HSTRY(:,:,valTIMESTEP + 1) = matCOEFF;
        
        if flgGIF == true
            hFig1 = fcnGIF(valTIMESTEP, valNELE, matDVE, matVLST, matCENTER, matELST, matDVECT, matPLEX, matCOEFF, matUINF, matROTANG, ...
                valWNELE, matWDVE, matWVLST, matWCENTER, matWELST, matWDVECT, vecWDVESURFACE, valWSIZE, matWVGRID, valGIFNUM);
        end
        
        %% Calculating Forces
        [matF_FS(:,:,valTIMESTEP), matF_IF(:,:,valTIMESTEP), matF_ID(:,:,valTIMESTEP), matINTCIRC(:,valTIMESTEP)] = ...
            fcnDVEFORCES2(matVLST, matELST, matROTANG, ...
            matCOEFF, vecTEDVE, valDENSITY, valNELE, matSPANDIR, vecTE, matPLEX, matWCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG, ...
            matWVLST, vecWLEDVE, matWELST, valWSIZE, matWDVE, vecWDVEFLIP, matWEIDX, vecWEMU, ...
            vecWVMU, matCENTER, matKINCON_P, vecKINCON_DVE, matUINF_KK, matDVECT, valZTOL, vecWDVESDFLIP, vecDVESDFLIP);
        
        [CL(valTIMESTEP), CDi(valTIMESTEP), e(valTIMESTEP), ~, ~, ~, ~, ~, CT(valTIMESTEP,:)] = fcnFORCES(valTIMESTEP, valWINGS, valROTORS, ...
            matF_FS, matF_IF, matF_ID, matLIFT_DIR, matSIDE_DIR, matDRAG_DIR, valDENSITY, valAREA, valSPAN, valUINF, vecROTORRPM, ...
            vecROTORDIAM, matROTORAXIS, vecTSITER, vecRSQUARED, vecDVEROTOR, vecDVEWING);
        
        matSPANDIR_ALL(:,:,valTIMESTEP) = matSPANDIR;
        matUINF_ALL(:,:,valTIMESTEP) = matUINF;
        
    end
    
    vecTSTIME(valTIMESTEP) = toc;
end

%% Apparent mass
% skip = 1;
% [CT_U, CL_U, matDGAMMADT, matF_NC] = fcnDGAMMADT(skip, valNELE, strATYPE, valDELTIME, matSPANDIR_ALL, matUINF_ALL, matF_FS, matF_IF, matF_ID, matINTCIRC, vecLIFT_DIR, vecSIDE_DIR, vecDRAG_DIR, valDENSITY, valAREA, valSPAN, valDIAM, valRPM);

%% Save environment
save(['run_',num2str(now),'.mat']);


