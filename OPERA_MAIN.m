clear

% profile('off')
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
% strFILE = 'inputs/ellipse.dat';
% strFILE = 'inputs/kussner.dat'
strFILE = 'inputs/TMotor_coarse2.dat'

[matPOINTS, strATYPE, vecSYM, flagRELAX, valMAXTIME, valDELTIME, valALPHA, ...
    valBETA, matTEPOINTS, matLEPOINTS, vecULS, valAREA, valSPAN, valDENSITY, vecDVESYM, valDIAM, ...
    valCOLL, valRPM, valJ, vecDVESURFACE, vecDVEFLIP, valBLADES, valGUSTAMP, valGUSTL, valGUSTMODE, valGUSTSTART, ...
    strHTYPE, matHINGELOC, matHINGEDIR, matHINGEANG, vecM] = fcnOPREAD(strFILE);

[TR, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matPLEX, matDVECT, matVATT, ~, matCENTER, matROTANG, ~, vecDVEAREA, matSPANDIR, vecDVECHORD]...
    = fcnTRIANG(matPOINTS, vecDVEFLIP);

% valJ = J(jj)

flagGIF = 1;
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

% Calculating matUINF for rotor or wing
if strcmpi(strATYPE{1}, 'ROTOR') || strcmpi(strATYPE{1}, 'PROPELLER')
    vecHUB = [0 0 0];
    vecROTORRADPS = valRPM.*2.*pi./60;
    
    if strcmpi(strATYPE{1}, 'ROTOR')
        translation = valJ.*(valRPM.*(pi/30)).*(valDIAM/2).*fcnUINFWING(valALPHA, 0);
    elseif strcmpi(strATYPE{1}, 'PROPELLER')
        translation = valJ.*(valRPM/60).*(valDIAM).*fcnUINFWING(valALPHA, 0);
    end
    matUINF = cross(repmat([0,0,-vecROTORRADPS],length(matCENTER(:,1)),1),matCENTER) - translation;
    matVUINF = cross(repmat([0,0,-vecROTORRADPS],length(matVLST(:,1)),1),matVLST) - translation;
else
    matUINF = repmat(fcnUINFWING(valALPHA, 0), valNELE, 1).*valUINF;
    matVUINF = repmat(fcnUINFWING(valALPHA, 0), size(matVLST,1), 1).*valUINF;
    matUINF_OR = matUINF;
    matAINF = matUINF.*0;
end

% Creating GIF file
if flagGIF == 1
    if exist('GIF','dir') ~= 7
        mkdir('GIF');
    end
    clear dir
    valGIFNUM = size(dir('GIF'),1) - 1;
end

% hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], matUINF, matROTANG, [], 'opengl');

%% D-Matrix Creation
matKINCON_DVE = repmat([1:valNELE]',4,1);
matKINCON_P = [matCENTER; 0.8.*(matVLST(matDVE(:,1),:) - matCENTER) + matCENTER; 0.8.*(matVLST(matDVE(:,2),:) - matCENTER) + matCENTER; 0.8.*(matVLST(matDVE(:,3),:) - matCENTER) + matCENTER];

% Points where flow tan-gency is enforced
matD = fcnDWING9(strATYPE, matEATT, matPLEX, valNELE, matELST, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecLEDVE, vecTEDVE, matROTANG, matSPANDIR, matKINCON_P, matKINCON_DVE, vecSYM, vecSYMDVE, vecDVESYM);
valDLEN = length(matD);

%% Preparing to timestep
% Initializing
matWELST = []; matWDVE = []; valWNELE = []; matWEATT = []; matWEIDX = []; matWPLEX = []; matWDVECT = [];
matWALIGN = []; matWCENTER = []; matWAKEGEOM = []; matWCOEFF = []; matWVLST = []; matWROTANG = [];
matWPLANE = []; vecWVMU = []; vecWEMU = []; matWVGRID = []; matWEGRID = []; vecWDVEFLIP = logical([]); vecWLEDVE = []; vecWTEDVE = [];
vecWDVECIRC = []; CL = nan(valMAXTIME,1); CDi = nan(valMAXTIME,1); CT = nan(valMAXTIME,1);
e = nan(valMAXTIME,1); gust_vel_old = matCENTER.*0; vecWDVESYM = logical([]); vecWSYMDVE = []; vecWSYM = []; valWSIZE = length(vecTE); vecWDVESURFACE = [];
vecDGAMMA_DT = zeros(size(vecDVESURFACE)); valPRESTEPS = 0; strWAKE_TYPE = strATYPE{3}; matINTCIRC = nan(1, length(vecTEDVE));
matLIFTFREE = []; matSIDEFREE = []; matLIFTIND = []; matSIDEIND = []; matDRAGIND = [];
matDVEDRAG_DIR = []; matDVELIFT_DIR = []; matDVESIDE_DIR = []; vecTSITER = nan(valMAXTIME,1); matDGAMMADT = []; matDVENC = []; vecDGAMMA_DETA = [];
vecRSQUARED = nan(valMAXTIME,2);

% Building wing resultant
vecR = fcnRWING(valDLEN, 0, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matKINCON_P, matKINCON_DVE, matDVECT, []);
boolKINCON = false(size(matD,1),1);
boolKINCON(end-(valNELE*4-1):end) = true;

% Solving for wing coefficients
matCOEFF = fcnSOLVED(matD, vecR, valNELE);
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
        [matVLST, matCENTER, matNEWWAKE, matUINF, matVUINF, matPLEX, matDVECT, matROTANG, vecHUB, matSPANDIR, matKINCON_P] = fcnMOVEROTOR(strATYPE{1}, valRPM, valJ, valDIAM, valALPHA, valDELTIME, matVLST, matELST, vecTE, matDVE, matCENTER, vecDVEFLIP, vecHUB, matSPANDIR, matKINCON_P);
    else
        [matVLST, matCENTER, matNEWWAKE, matKINCON_P] = fcnMOVEWING(matUINF_OR, valDELTIME, matVLST, matCENTER, matELST, vecTE, matKINCON_P);
    end
    
    if ~isempty(valGUSTSTART) && valGUSTSTART >= 0 && valTIMESTEP > valPRESTEPS
        tmatUINF = matUINF;
        [matUINF, gust_vel_old] = fcnGUSTWING(matUINF, valGUSTAMP, valGUSTL, valGUSTMODE, valDELTIME, valUINF, valGUSTSTART, matCENTER, gust_vel_old, valPRESTEPS*valDELTIME*10);
        matAINF = (matUINF - tmatUINF)./valDELTIME;
    end
    
    % Generating new wake elements
    if any(vecTE)
        [matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWPLEX, matWDVECT, matWCENTER, matWCOEFF, matWROTANG, vecWLE, vecWTE, ...
            vecWLEDVE, vecWTEDVE, vecWDVECIRC, vecWSYM, vecWSYMDVE, vecWDVESYM, matWVGRID, vecWOTE, ...
            vecWOTEDVE, matWEGRID, matWE2GRID, vecWVMU, vecWEMU, vecWDVEFLIP, vecWDVESURFACE] = fcnCREATEWAKE2(valTIMESTEP, strATYPE, matNEWWAKE, matCOEFF, valWSIZE, ...
            vecTEDVE, matCENTER, matROTANG, matWCOEFF, matWPLEX, vecWDVECIRC, vecWSYMDVE, vecSYMDVE, vecWDVESYM, vecDVESYM, ...
            vecWSYM, matWVGRID, matWVLST, matWELST, matWEGRID, vecWVMU, vecWEMU, vecWDVEFLIP, matWCENTER, matWROTANG, matWDVECT, ...
            matWDVE, matWEIDX, vecWLEDVE, matWEATT, vecDVESURFACE, vecWDVESURFACE);
        
        if valTIMESTEP > valPRESTEPS
            
            count = 0;
            delt = 1;
            delt_old = 5;
            while ~isnan(delt) && delt >= 0.001
                int_circ1 = fcnINTCIRC2(matPLEX, matCOEFF, matDVEGRID);

                [vecR, wind] = fcnRWING(valDLEN, valTIMESTEP, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matKINCON_P, matKINCON_DVE, matDVECT, vecWDVESYM);
                matCOEFF = fcnSOLVED(matD, vecR, valNELE);

                % Update wake coefficients
                [vecWVMU, vecWEMU] = fcnWAKEMU(strATYPE, vecWLE, matWVGRID, matWEGRID, matWE2GRID, vecWVMU, vecWEMU, matWELST, matWVLST, vecTEDVE, matCOEFF, matCENTER, matROTANG, vecWOTE);
                matWCOEFF = fcnADJCOEFF(vecWVMU, vecWEMU, matWVLST, matWCENTER, matWROTANG, matWDVE, matWCOEFF, matWELST, matWEIDX, valWNELE);
                
                int_circ2 = fcnINTCIRC2(matPLEX, matCOEFF, matDVEGRID);
                delt = abs((int_circ2 - int_circ1)./int_circ1);
                
                count = count + 1;
                
                if delt_old < delt
                    disp(['          Diverged. Error: ', num2str(delt), ' Count = ', num2str(count)]);
                    break;
                end
                
                delt_old = delt;
            end    
            vecTSITER(valTIMESTEP) = count;
            [vecRSQUARED(valTIMESTEP,1), vecRSQUARED(valTIMESTEP,2)] = fcnRSQUARED(matD, vecR, matCOEFF, boolKINCON, matVLST, matVATT, matCENTER, matROTANG, matELST, matEATT);
            
            % Relaxing Wake
            if flagRELAX == 1
                [matWELST, matWVLST, matWDVE, valWNELE, matWEIDX, matWPLEX, matWDVECT, matWCENTER, matWROTANG, matWVGRID] = ...
                    fcnRELAX7(flagHVRMOD, valDELTIME, valTIMESTEP, valRPM, valNELE, matCOEFF, matPLEX, valWNELE, matWCOEFF, matWDVE, matWVLST, matWPLEX, valWSIZE, ...
                    matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST, matWEIDX, vecDVESYM, vecWDVESYM, vecWSYM, matWVGRID, vecWDVEFLIP, valPRESTEPS, matWE2GRID);
                
                % Update all wake coefficients
                matWCOEFF = fcnADJCOEFF(vecWVMU, vecWEMU, matWVLST, matWCENTER, matWROTANG, matWDVE, matWCOEFF, matWELST, matWEIDX, valWNELE);
            end
            
            % Updating coefficient convergence history
            matCOEFF_HSTRY(:,:,valTIMESTEP + 1) = matCOEFF;
            
            if flagGIF == 1
                hFig1 = fcnGIF(valTIMESTEP - valPRESTEPS, valNELE, matDVE, matVLST, matCENTER, matELST, matDVECT, matPLEX, matCOEFF, matUINF, matROTANG, ...
                    valWNELE, matWDVE, matWVLST, matWCENTER, matWELST, matWDVECT, vecWDVESURFACE, valWSIZE, 0, matWVGRID, valGIFNUM);
            end
            
            %% Calculating Forces
            [vecDVELIFT, vecDVEDRAG, vecDVESIDE, matDVEDRAG_DIR, matDVELIFT_DIR, matDVESIDE_DIR, vecDGAMMA_DT, vecDGAMMA_DETA, matINTCIRC, matLIFTFREE, matSIDEFREE, matLIFTIND, matSIDEIND, matDRAGIND] = fcnDVEFORCES(strATYPE{3},valTIMESTEP, strWAKE_TYPE, matVLST, matELST, matROTANG, ...
                matUINF, matCOEFF, vecTEDVE, valDENSITY, valNELE, matSPANDIR, vecTE, matPLEX, matWCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG, ...
                matVUINF, matWVLST, vecWLE, vecWLEDVE, matWELST, valAREA, valSPAN, vecDVESYM, vecWDVESYM, vecDGAMMA_DT, vecDGAMMA_DETA, valWSIZE, matWDVE, vecWDVEFLIP, matWEIDX, vecWEMU, vecWVMU, matINTCIRC, matLIFTFREE, matSIDEFREE, matLIFTIND, matSIDEIND, matDRAGIND, ...
                matDVEDRAG_DIR, matDVELIFT_DIR, matDVESIDE_DIR, matDVEGRID, vecDVEAREA, vecCHORD);
            
            if strcmpi(strATYPE{1}, 'WING')
                [CL(valTIMESTEP), CDi(valTIMESTEP), CY(valTIMESTEP), e(valTIMESTEP)] = fcnWFORCES(valTIMESTEP, vecDVELIFT, vecDVEDRAG, vecDVESIDE, valDENSITY, valAREA, valSPAN, vecTSITER, vecRSQUARED);
            else
                [CT(valTIMESTEP), vecDVETHRUST] = fcnPFORCES(strATYPE{1}, valTIMESTEP, vecDVELIFT, vecDVEDRAG, vecDVESIDE, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR, valDENSITY, valDIAM, valRPM, vecTSITER, vecRSQUARED);
            end
            
            %% Unsteady
            matSPANDIR_ALL(:,:,valTIMESTEP) = matSPANDIR;
            matUINF_ALL(:,:,valTIMESTEP) = matUINF;
            
        end
    end
    vecTSTIME(valTIMESTEP) = toc;
end

%%
% Adding in apparent mass
skip = 1;
[CT_U, CL_U, matDGAMMADT, matDVENC] = fcnDGAMMADT(skip, valDELTIME, strATYPE, matINTCIRC, valDENSITY, valRPM, valDIAM, valAREA, valUINF, matLIFTFREE, matLIFTIND, matDRAGIND, matSIDEFREE, matSIDEIND, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR, matSPANDIR_ALL, matUINF_ALL, vecTE, vecTEDVE);

save(['tmotor_run_',num2str(now),'.mat']);

figure(20);
plot(CT, '-k')
hold on
plot(CT_U, '--b');
hold off
grid minor
box on
axis tight

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
% plot(s_t - 1/(c/2), CLOG2D, 'rs');
% hold off
% xlim([0 8])
