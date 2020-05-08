% clear

% profile clear
% profile('-memory','on');

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
% strFILE = 'inputs/test.dat'
% strFILE = 'inputs/ellipse.dat';
% strFILE = 'inputs/kussner.dat'
strFILE = 'inputs/TMotor_coarse.dat'

[matPOINTS, strATYPE, vecSYM, flagRELAX, valMAXTIME, valDELTIME, valALPHA, ...
    valBETA, matTEPOINTS, matLEPOINTS, vecULS, valAREA, valSPAN, valDENSITY, vecDVESYM, valDIAM, ...
    valCOLL, valRPM, valJ, vecDVESURFACE, vecDVEFLIP, valBLADES, valGUSTAMP, valGUSTL, valGUSTMODE, valGUSTSTART, ...
    strHTYPE, matHINGELOC, matHINGEDIR, matHINGEANG, vecM] = fcnOPREAD(strFILE);

[TR, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matPLEX, matDVECT, matVATT, ~, matCENTER, matROTANG, ~, vecDVEAREA, matSPANDIR, vecDVECHORD]...
    = fcnTRIANG(matPOINTS, vecDVEFLIP);

% valALPHA = 15
% valJ = 0.2889
valALPHA = SEQ_ALPHA(jjj)
valJ = SEQ_JJ(jjj)

flagRELAX = 1;
flagWALL = false;

valMAXTIME = 320;
valDELTIME = 0.00025;

valSTARTFORCES = 0;

flagGIF = 0;

flagHVRMOD = false;
valUINF = 1;

%% Preliminary Geometry Stuff
% Duplicating blades for rotors
if valBLADES > 1
    [valNELE, matVLST, matELST, matDVE, matCENTER, matPLEX, vecDVEAREA, ...
        matEATT, matEIDX, vecDVESURFACE, vecDVEFLIP, vecDVESYM, matDVECT, ...
        matROTANG, matVATT, matTEPOINTS, matLEPOINTS, matSPANDIR, matHINGELOC, matHINGEDIR] = fcnDUPBLADE(valBLADES, valNELE, matVLST, matELST, matDVE, ...
        matCENTER, matPLEX, vecDVEAREA, matEATT, matEIDX, vecDVESURFACE, vecDVEFLIP, vecDVESYM, matDVECT, matROTANG, matVATT, matTEPOINTS, matLEPOINTS, matSPANDIR, matHINGELOC, matHINGEDIR);
end

% Generating leading and trailing edge information
[vecLE, vecLEDVE, vecTE, vecTEDVE, vecSYM, vecSYMDVE, matELST, matDVEGRID, vecCHORD, vecSPAN] = fcnLETEGEN(strATYPE, matVLST, matELST, matEATT, matLEPOINTS, matTEPOINTS, vecSYM, matEIDX, vecM);

if ~isempty(strHTYPE)
    % Setting blade angles for propellers
    [matHINGELOC, matHINGEDIR, matVLST, matCENTER, matPLEX, matDVECT, matROTANG, matSPANDIR] = fcnSETHANG(matHINGEANG, matHINGELOC, matHINGEDIR, matVLST, matDVE, vecDVESURFACE, valBLADES, vecDVEFLIP, matSPANDIR);
end

% Kinematic condition locations
tmp = find(vecDVEAREA > .5e-4);
matKINCON_DVE = [1:valNELE repmat(tmp',1,3)]';
dist = 0.8;
matKINCON_P = [matCENTER; dist.*(matVLST(matDVE(tmp,1),:) - matCENTER(tmp,:)) + matCENTER(tmp,:); dist.*(matVLST(matDVE(tmp,2),:) - matCENTER(tmp,:)) + matCENTER(tmp,:); dist.*(matVLST(matDVE(tmp,3),:) - matCENTER(tmp,:)) + matCENTER(tmp,:)];

% Calculating matUINF for rotor or wing
if strcmpi(strATYPE{1}, 'ROTOR') || strcmpi(strATYPE{1}, 'PROPELLER')
    valAZPERREV = (1/(valRPM/60))/valDELTIME;
    if mod(valAZPERREV, 1) ~= 0
        valAZPERREV = ceil((1/(valRPM/60))/valDELTIME);
        fprintf('\nTimestep size corrected from %g to %g (%d azimuth stations per revolution)\n\n', valDELTIME, (1/(valRPM/60))/valAZPERREV, valAZPERREV)
        valDELTIME = (1/(valRPM/60))/valAZPERREV;
    end
    
    vecHUB = [0 0 0];
    vecROTORRADPS = valRPM.*2.*pi./60;
    
    if strcmpi(strATYPE{1}, 'ROTOR')
        translation = valJ.*(valRPM.*(pi/30)).*(valDIAM/2).*fcnUINFWING(valALPHA, 0);
    elseif strcmpi(strATYPE{1}, 'PROPELLER')
        translation = valJ.*(valRPM/60).*(valDIAM).*fcnUINFWING(valALPHA, 0);
    end
    matUINF = cross(repmat([0,0,-vecROTORRADPS],length(matCENTER(:,1)),1),matCENTER) - translation;
    matVUINF = cross(repmat([0,0,-vecROTORRADPS],length(matVLST(:,1)),1),matVLST) - translation;
    matUINF_KK = cross(repmat([0,0,-vecROTORRADPS],length(matKINCON_P(:,1)),1),matKINCON_P) - translation;
    vecLIFT_DIR = [];
    vecSIDE_DIR = [];
    vecDRAG_DIR = [];
    
    if flagWALL == true
        valWOFF = 0.45; % wall offset, meters from the hub to the wall (in the X-Z plane)
        vecWOFF = [sind(valALPHA) 0 -cosd(valALPHA)];
        vecWNORM = [-sind(valALPHA) 0 cosd(valALPHA)];
    else
        valWOFF = nan;
        vecWOFF = [nan nan nan];
        vecWNORM = [nan nan nan];
    end
    
else
    matUINF = repmat(fcnUINFWING(valALPHA, 0), valNELE, 1).*valUINF;
    matVUINF = repmat(fcnUINFWING(valALPHA, 0), size(matVLST,1), 1).*valUINF;
    matUINF_KK = repmat(fcnUINFWING(valALPHA, 0), size(matKINCON_P,1), 1).*valUINF;
    valAZPERREV = nan;
    
    vecDRAG_DIR = matUINF(1,:)./sqrt(sum(matUINF(1,:).^2,2));
    vecLIFT_DIR = cross(vecDRAG_DIR, matSPANDIR(1,:), 2);
    vecLIFT_DIR = vecLIFT_DIR./sqrt(sum(vecLIFT_DIR.^2,2));
    vecSIDE_DIR = cross(vecLIFT_DIR, vecDRAG_DIR, 2);
    vecSIDE_DIR = vecSIDE_DIR./sqrt(sum(vecSIDE_DIR.^2,2));
end

% Creating GIF file
if flagGIF == 1
    if exist('GIF','dir') ~= 7
        mkdir('GIF');
    end
    clear dir
    valGIFNUM = size(dir('GIF'),1) - 1;
end

%% D-Matrix Creation
matD = fcnDWING9(strATYPE, matEATT, matPLEX, valNELE, matELST, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecLEDVE, vecTEDVE, matROTANG, matSPANDIR, matKINCON_P, matKINCON_DVE, vecSYM, vecSYMDVE, vecDVESYM);
valDLEN = length(matD);

%% Preparing to timestep
valZTOL = 1e-4; % Used to offset field points from surface sheets when necessary

% Initializing
matWELST = []; matWDVE = []; valWNELE = []; matWEATT = []; matWEIDX = []; matWPLEX = []; matWDVECT = [];
matWALIGN = []; matWCENTER = []; matWAKEGEOM = []; matWCOEFF = []; matWVLST = []; matWROTANG = [];
matWPLANE = []; vecWVMU = []; vecWEMU = []; matWVGRID = []; matWEGRID = []; vecWDVEFLIP = logical([]); vecWLEDVE = []; vecWTEDVE = [];
vecWDVECIRC = []; CL = nan(valMAXTIME,1); CDi = nan(valMAXTIME,1); CT = nan(valMAXTIME,1);
e = nan(valMAXTIME,1); gust_vel_old = matKINCON_P.*0; vecWDVESYM = logical([]); vecWSYMDVE = []; vecWSYM = []; valWSIZE = length(vecTE); vecWDVESURFACE = [];
vecDGAMMA_DT = zeros(size(vecDVESURFACE)); valPRESTEPS = 0; strWAKE_TYPE = strATYPE{3}; matINTCIRC = nan(valNELE,valMAXTIME);
matLIFTFREE = []; matSIDEFREE = []; matLIFTIND = []; matSIDEIND = []; matDRAGIND = [];
matDVEDRAG_DIR = []; matDVELIFT_DIR = []; matDVESIDE_DIR = []; vecTSITER = nan(valMAXTIME,1); matDGAMMADT = []; matDVENC = []; vecDGAMMA_DETA = [];
vecRSQUARED = nan(valMAXTIME,2); matWDVEGRID = [];

% Building wing resultant
vecR = fcnRWING(strATYPE, valZTOL, valDLEN, 0, matUINF_KK, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matKINCON_P, matKINCON_DVE, matDVECT, vecWDVESYM, matWDVECT, []);
boolKINCON = false(size(matD,1),1);
boolKINCON(end-(valNELE*4-1):end) = true;

% Solving for wing coefficients
matCOEFF = fcnSOLVED(matD, vecR, valNELE);
[vecVMU, vecEMU] = fcnVEMU(matVLST, matVATT, matCENTER, matROTANG, matCOEFF, matELST, matEATT, vecTE);
matCOEFF_HSTRY(:,:,1) = matCOEFF;

% Setting prestep conditions
if strcmpi(strATYPE{1}, 'WING')
    valPRESTEPS = 0;
    valDELTIME = valDELTIME*10;
    strWAKE_TYPE = 'STEADY';
elseif (strcmpi(strATYPE{1}, 'PROPELLER') && valJ == 0) || (strcmpi(strATYPE{1}, 'ROTOR') && valJ == 0)
    flagHVRMOD = true;
end

%% Timestep to solution
for valTIMESTEP = 1:valMAXTIME
    tic
    if strcmpi(strATYPE{1}, 'WING') && valTIMESTEP == valPRESTEPS + 1
        valDELTIME = valDELTIME/10;
        strWAKE_TYPE = strATYPE{3};
    end
    
    if (strcmpi(strATYPE{1}, 'ROTOR') || strcmpi(strATYPE{1}, 'PROPELLER'))
        [matVLST, matCENTER, matNEWWAKE, matUINF, matVUINF, matUINF_KK, matPLEX, matDVECT, matROTANG, vecHUB, matSPANDIR, matKINCON_P] = fcnMOVEROTOR(strATYPE{1}, valRPM, valJ, valDIAM, valALPHA, valDELTIME, matVLST, matELST, vecTE, matDVE, matCENTER, vecDVEFLIP, vecHUB, matSPANDIR, matKINCON_P);
    else
        [matVLST, matCENTER, matNEWWAKE, matKINCON_P] = fcnMOVEWING(matUINF, valDELTIME, matVLST, matCENTER, matELST, vecTE, matKINCON_P);
    end
    
    if ~isempty(valGUSTSTART) && valGUSTSTART >= 0 && valTIMESTEP > valPRESTEPS && strcmpi(strATYPE{1}, 'WING')
        [matUINF_KK, gust_vel_old] = fcnGUSTWING(matUINF_KK, valGUSTAMP, valGUSTL, valGUSTMODE, valDELTIME, valUINF, valGUSTSTART, matKINCON_P, gust_vel_old, valPRESTEPS*valDELTIME*10);
    end
    
    % Generating new wake elements
    if any(vecTE)
        [matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWPLEX, matWDVECT, matWCENTER, matWCOEFF, matWROTANG, vecWLE, vecWTE, ...
            vecWLEDVE, vecWTEDVE, vecWDVECIRC, vecWSYM, vecWSYMDVE, vecWDVESYM, matWVGRID, vecWOTE, ...
            vecWOTEDVE, matWEGRID, matWE2GRID, vecWVMU, vecWEMU, vecWDVEFLIP, vecWDVESURFACE, matWDVEGRID] = fcnCREATEWAKE2(valTIMESTEP, strATYPE, matNEWWAKE, matCOEFF, valWSIZE, ...
            vecTEDVE, matCENTER, matROTANG, matWCOEFF, matWPLEX, vecWDVECIRC, vecWSYMDVE, vecSYMDVE, vecWDVESYM, vecDVESYM, ...
            vecWSYM, matWVGRID, matWVLST, matWELST, matWEGRID, vecWVMU, vecWEMU, vecWDVEFLIP, matWCENTER, matWROTANG, matWDVECT, ...
            matWDVE, matWEIDX, vecWLEDVE, matWEATT, vecDVESURFACE, vecWDVESURFACE, valAZPERREV, matWDVEGRID);
        
        if valTIMESTEP > valPRESTEPS
            
            [vecR, matCOEFF, vecVMU, vecEMU, matWCOEFF, vecWVMU, vecWEMU, vecTSITER, vecRSQUARED] = ...
                fcnTSITER(matPLEX, matCOEFF, matDVEGRID, valDLEN, valTIMESTEP, matUINF_KK, valWNELE, ...
                matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matKINCON_P, matKINCON_DVE, ...
                matDVECT, vecWDVESYM, matD, vecR, valNELE, matVLST, matVATT, matCENTER, matROTANG, ...
                matDVE, matELST, matEIDX, strATYPE, vecWLE, matWVGRID, matWEGRID, matWE2GRID, ...
                vecWVMU, vecWEMU, matWELST, matWVLST, vecTEDVE, vecWOTE, matWDVE, matWEIDX, ...
                matEATT, boolKINCON, vecTE, vecVMU, vecEMU, valAZPERREV, matWDVECT, valZTOL, vecRSQUARED);
            
            % Relaxing Wake
            if flagRELAX == 1
%                 hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, [], 'opengl');
%                 hFig1 = fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, valWSIZE, valPRESTEPS, matWVGRID, vecWDVESURFACE);
%                 view([-17 29]);
                [matWELST, matWVLST, matWPLEX, matWDVECT, matWCENTER, matWROTANG] = ...
                    fcnRELAX9(valTIMESTEP, valDELTIME, valNELE, matCOEFF, matPLEX, valWNELE, matWCOEFF, matWDVE, matWVLST, matWPLEX, valWSIZE, ...
                    matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST, matWEIDX, vecDVESYM, vecWDVESYM, matWVGRID, vecWDVEFLIP, ...
                    valPRESTEPS, matWDVEGRID, valWOFF, vecWOFF, vecWNORM, vecHUB, matVLST, matDVE, vecDVEFLIP);
                
                matWCOEFF = fcnADJCOEFF(vecWVMU, vecWEMU, matWVLST, matWCENTER, matWROTANG, matWDVE, matWCOEFF, matWELST, matWEIDX, valWNELE);
                
                [vecR, matCOEFF, vecVMU, vecEMU, matWCOEFF, vecWVMU, vecWEMU, vecTSITER, vecRSQUARED] = ...
                    fcnTSITER(matPLEX, matCOEFF, matDVEGRID, valDLEN, valTIMESTEP, matUINF_KK, valWNELE, ...
                    matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matKINCON_P, matKINCON_DVE, ...
                    matDVECT, vecWDVESYM, matD, vecR, valNELE, matVLST, matVATT, matCENTER, matROTANG, ...
                    matDVE, matELST, matEIDX, strATYPE, vecWLE, matWVGRID, matWEGRID, matWE2GRID, ...
                    vecWVMU, vecWEMU, matWELST, matWVLST, vecTEDVE, vecWOTE, matWDVE, matWEIDX, ...
                    matEATT, boolKINCON, vecTE, vecVMU, vecEMU, valAZPERREV, matWDVECT, valZTOL, vecRSQUARED);
            end
            
            % Updating coefficient convergence history
            matCOEFF_HSTRY(:,:,valTIMESTEP + 1) = matCOEFF;
            
            if flagGIF == 1
                hFig1 = fcnGIF(valTIMESTEP - valPRESTEPS, valNELE, matDVE, matVLST, matCENTER, matELST, matDVECT, matPLEX, matCOEFF, matUINF, matROTANG, ...
                    valWNELE, matWDVE, matWVLST, matWCENTER, matWELST, matWDVECT, vecWDVESURFACE, valWSIZE, 0, matWVGRID, valGIFNUM);
            end
            
            %% Calculating Forces
            %             hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, [], 'opengl');
            if valTIMESTEP > valSTARTFORCES
                [matF_FS(:,:,valTIMESTEP), matF_IF(:,:,valTIMESTEP), matF_ID(:,:,valTIMESTEP), matINTCIRC(:,valTIMESTEP)] = ...
                    fcnDVEFORCES2(valTIMESTEP, matVLST, matELST, matROTANG, ...
                    matCOEFF, vecTEDVE, valDENSITY, valNELE, matSPANDIR, vecTE, matPLEX, matWCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG, ...
                    matWVLST, vecWLEDVE, matWELST, vecDVESYM, vecWDVESYM, valWSIZE, matWDVE, vecWDVEFLIP, matWEIDX, vecWEMU, ...
                    vecWVMU, matCENTER, matKINCON_P, matKINCON_DVE, matUINF_KK, vecWLE, matWDVECT, matDVECT, valZTOL);
            else
                matF_FS(:,:,valTIMESTEP) = nan(valNELE,3);
                matF_IF(:,:,valTIMESTEP) = nan(valNELE,3);
                matF_ID(:,:,valTIMESTEP) = nan(valNELE,3);
            end
            
            if strcmpi(strATYPE{1}, 'WING')
                [CL(valTIMESTEP), CDi(valTIMESTEP), e(valTIMESTEP), ~,~,~,~,~] = fcnWFORCES(true, valTIMESTEP, matF_FS, matF_IF, matF_ID, vecLIFT_DIR, vecSIDE_DIR, vecDRAG_DIR, valDENSITY, valAREA, valSPAN, vecTSITER, vecRSQUARED);
            else
                [CT(valTIMESTEP), vecDVETHRUST] = fcnPFORCES(true, strATYPE{1}, valTIMESTEP, matF_FS, matF_IF, matF_ID, valDENSITY, valDIAM, valRPM, vecTSITER, vecRSQUARED);
            end
            
            
            %% Unsteady
            matSPANDIR_ALL(:,:,valTIMESTEP) = matSPANDIR;
            matUINF_ALL(:,:,valTIMESTEP) = matUINF;
            
        end
    end
    vecTSTIME(valTIMESTEP) = toc;
    
    %     deg_per_ts = valRPM.*(pi/30).*(180/pi).*valDELTIME;
    %     pos = valTIMESTEP.*deg_per_ts + 90;
    %     pos = mod(pos,360);
    %
    %     hFig1 = fcnGIF(valTIMESTEP, valNELE, matDVE, matVLST, matCENTER, matELST, matDVECT, matPLEX, matCOEFF, matUINF, matROTANG, ...
    %         valWNELE, matWDVE, matWVLST, matWCENTER, matWELST, matWDVECT, vecWDVESURFACE, valWSIZE, 0, matWVGRID, valGIFNUM, pos);
    
    save('tmp.mat');
end

%% Apparent mass
skip = 1;
[CT_U, CL_U, matDGAMMADT] = fcnDGAMMADT(skip, valNELE, strATYPE, valDELTIME, matSPANDIR_ALL, matUINF_ALL, matF_FS, matF_IF, matF_ID, matINTCIRC, vecLIFT_DIR, vecSIDE_DIR, vecDRAG_DIR, valDENSITY, valAREA, valSPAN, valDIAM, valRPM);

%%


save(['tmotor_run_',num2str(now),'.mat']);

% figure(20);
% plot(CT, '-k')
% hold on
% plot(CT_U, '--b');
% hold off
% grid minor
% box on
% axis tight

% profile report

% load('Kussner.mat')
% hFig25 = figure(25);
% clf(25);
% plot(((s(1:40))), c_l(1:40), '-k');
% box on
% axis tight
% grid minor
% hold on
% 
% c = 1
% s_t = ([1:skip:valMAXTIME].*valDELTIME.*2)/c;
% valAR = (valSPAN.^2)./valAREA;
% CLOG2D = CL_U.*((valAR + 2)/valAR);
% plot(s_t - (valGUSTSTART*valDELTIME)/(c/2), CL(1:skip:end).*((valAR + 2)/valAR), '-.md')
% plot(s_t - (valGUSTSTART*valDELTIME)/(c/2), CLOG2D, 'rs');
% hold off
% xlim([0 8])




