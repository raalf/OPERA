clear

%profile('off')
%profile('-memory','on');

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
% strFILE = 'inputs/goland_wing.dat'
% strFILE = 'inputs/kussner.dat'
% strFILE = 'inputs/box_wing.dat'
% strFILE = 'inputs/TMotor.dat'
strFILE = 'inputs/TMotor_coarse.dat'
% strFILE = 'inputs/TMotor_nocamber.dat'
% strFILE = 'inputs/Leishman_Rotor.dat'
% strFILE = 'inputs/Caradonna_Rotor.dat'
% strFILE = 'Stuff/Columbia/Columbia_234_Rotor.dat'
% strFILE = 'inputs/test.dat'
% strFILE = 'inputs/MAE_11x7.dat'

[matPOINTS, strATYPE, vecSYM, flagRELAX, valMAXTIME, valDELTIME, valALPHA, ...
    valBETA, matTEPOINTS, matLEPOINTS, vecULS, valAREA, valSPAN, valDENSITY, vecDVESYM, valDIAM, ...
    valCOLL, valRPM, valJ, vecDVESURFACE, vecDVEFLIP, valBLADES, valGUSTAMP, valGUSTL, valGUSTMODE, valGUSTSTART, ...
    strHTYPE, matHINGELOC, matHINGEDIR, matHINGEANG, vecM] = fcnOPREAD(strFILE);

[TR, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matPLEX, matDVECT, matVATT, ~, matCENTER, matROTANG, ~, vecDVEAREA, matSPANDIR, vecDVECHORD]...
    = fcnTRIANG(matPOINTS, vecDVEFLIP);

flagGIF = 0;
flagHVRMOD = false;
valUINF = 1;

% flagRELAX = 0
% valJ = J(jj)

% % Hinge
% strHTYPE = 'FPL';
% matHINGELOC = cat(3, [0 0.2032 0], [0 0.2032 0], [0 0.7493 0]);
% matHINGEDIR = cat(3, [1 0 0], [0 1 0], [0 0 1]);
% matHINGEANG = repmat([FLAP PITCH 0], 3, 1);

% % Hinge
% strHTYPE = 'FPL';
% matHINGELOC = cat(3, [0 0 0], [0 0 0], [0 0 0]);
% matHINGEDIR = cat(3, [1 0 0], [0 1 0], [0 0 1]);
% % matHINGEANG = repmat([0 -2.0 0], 3, 1);
% matHINGEANG = repmat([0 0 0], 3, 1);

%% Preliminary Geometry Stuff
% Duplicating blades for rotors
if valBLADES > 1
    [valNELE, matVLST, matELST, matDVE, matCENTER, matPLEX, vecDVEAREA, ...
        matEATT, matEIDX, vecDVESURFACE, vecDVEFLIP, vecDVESYM, matDVECT, ...
        matROTANG, matVATT, matTEPOINTS, matLEPOINTS, matSPANDIR, matHINGELOC, matHINGEDIR] = fcnDUPBLADE(valBLADES, valNELE, matVLST, matELST, matDVE, ...
        matCENTER, matPLEX, vecDVEAREA, matEATT, matEIDX, vecDVESURFACE, vecDVEFLIP, vecDVESYM, matDVECT, matROTANG, matVATT, matTEPOINTS, matLEPOINTS, matSPANDIR, matHINGELOC, matHINGEDIR);
end

% Generating leading and trailing edge information
[vecLE, vecLEDVE, vecTE, vecTEDVE, vecSYM, vecSYMDVE, matELST, matDVEGRID] = fcnLETEGEN(strATYPE, matVLST, matELST, matEATT, matLEPOINTS, matTEPOINTS, vecSYM, matEIDX, vecM);

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
matKINCON_DVE = [1:valNELE]';

% Points where flow tangency is enforced
matD = fcnDWING9(strATYPE, matEATT, matPLEX, valNELE, matELST, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecLEDVE, vecTEDVE, matROTANG, matSPANDIR, matCENTER, matKINCON_DVE, vecSYM, vecSYMDVE, vecDVESYM);
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
matDVEDRAG_DIR = []; matDVELIFT_DIR = []; matDVESIDE_DIR = [];

% Building wing resultant
vecR = fcnRWING(valDLEN, 0, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matCENTER, matKINCON_DVE, matDVECT, []);

% Solving for wing coefficients
matCOEFF = fcnSOLVED(matD, vecR, valNELE);
[vecVMU, vecEMU] = fcnVEMU(matVLST, matVATT, matCENTER, matROTANG, matCOEFF, matELST, matEATT, vecTE);
matCOEFF = fcnADJCOEFF(vecVMU, vecEMU, matVLST, matCENTER, matROTANG, matDVE, matCOEFF, matELST, matEIDX, valNELE);
matCOEFF_HSTRY(:,:,1) = matCOEFF;
vecDGAMMA_DETA = matCOEFF(:,2);

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
    if strcmpi(strATYPE{1}, 'WING') && valTIMESTEP == valPRESTEPS + 1
        valDELTIME = valDELTIME/10;
        strWAKE_TYPE = strATYPE{3};
    end
    
    if ~isempty(valGUSTSTART) && valGUSTSTART >= 0 && valTIMESTEP > valPRESTEPS
        tmatUINF = matUINF;
        [matUINF, gust_vel_old] = fcnGUSTWING(matUINF, valGUSTAMP, valGUSTL, valGUSTMODE, valDELTIME, valUINF, valGUSTSTART, matCENTER, gust_vel_old, valPRESTEPS*valDELTIME*10);
        matAINF = (matUINF - tmatUINF)./valDELTIME;
    end
    
    if (strcmpi(strATYPE{1}, 'ROTOR') || strcmpi(strATYPE{1}, 'PROPELLER'))
        [matVLST, matCENTER, matNEWWAKE, matUINF, matVUINF, matPLEX, matDVECT, matROTANG, vecHUB, matSPANDIR] = fcnMOVEROTOR(strATYPE{1}, valRPM, valJ, valDIAM, valALPHA, valDELTIME, matVLST, matELST, vecTE, matDVE, matCENTER, vecDVEFLIP, vecHUB, matSPANDIR);
    else
        [matVLST, matCENTER, matNEWWAKE] = fcnMOVEWING(matUINF_OR, valDELTIME, matVLST, matCENTER, matELST, vecTE);
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
            vecR = fcnRWING(valDLEN, valTIMESTEP, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matCENTER, matKINCON_DVE, matDVECT, vecWDVESYM);
            matCOEFF = fcnSOLVED(matD, vecR, valNELE);
            [vecVMU, vecEMU] = fcnVEMU(matVLST, matVATT, matCENTER, matROTANG, matCOEFF, matELST, matEATT, vecTE);
            matCOEFF = fcnADJCOEFF(vecVMU, vecEMU, matVLST, matCENTER, matROTANG, matDVE, matCOEFF, matELST, matEIDX, valNELE);
            
            % Update wake coefficients
            [vecWVMU, vecWEMU] = fcnWAKEMU(strATYPE, vecWLE, matWVGRID, matWEGRID, matWE2GRID, vecWVMU, vecWEMU, matWELST, matWVLST, vecTEDVE, matCOEFF, matCENTER, matROTANG, vecWOTE);
            matWCOEFF = fcnADJCOEFF(vecWVMU, vecWEMU, matWVLST, matWCENTER, matWROTANG, matWDVE, matWCOEFF, matWELST, matWEIDX, valWNELE);
            
            % Relaxing Wake
            if flagRELAX == 1
                [matWELST, matWVLST, matWDVE, valWNELE, matWEIDX, matWPLEX, matWDVECT, matWCENTER, matWROTANG, matWVGRID] = ...
                    fcnRELAX7(flagHVRMOD, valDELTIME, valTIMESTEP, valRPM, valNELE, matCOEFF, matPLEX, valWNELE, matWCOEFF, matWDVE, matWVLST, matWPLEX, valWSIZE, ...
                    matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST, matWEIDX, vecDVESYM, vecWDVESYM, vecWSYM, matWVGRID, vecWDVEFLIP, valPRESTEPS, matWE2GRID);
                
%                 [matWELST, matWVLST, matWDVE, valWNELE, matWEIDX, matWPLEX, matWDVECT, matWCENTER, matWROTANG, matWVGRID] = ...
%                     fcnRELAX8(flagHVRMOD, valDELTIME, valTIMESTEP, valRPM, valNELE, matCOEFF, matPLEX, valWNELE, matWCOEFF, matWDVE, matWVLST, matWPLEX, valWSIZE, ...
%                     matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST, matWEIDX, vecDVESYM, vecWDVESYM, vecWSYM, matWVGRID, vecWDVEFLIP, valPRESTEPS, matWE2GRID);
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
            %     hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, [], 'opengl');
            %     hFig1 = fcnPLOTWAKE(1, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, valWSIZE, valPRESTEPS, matWVGRID, vecWDVESURFACE);
            
            [vecDVELIFT, vecDVEDRAG, vecDVESIDE, matDVEDRAG_DIR, matDVELIFT_DIR, matDVESIDE_DIR, vecDGAMMA_DT, vecDGAMMA_DETA, matINTCIRC, matLIFTFREE, matSIDEFREE, matLIFTIND, matSIDEIND, matDRAGIND] = fcnDVEFORCES(valTIMESTEP, strWAKE_TYPE, matVLST, matELST, matROTANG, ...
                matUINF, matCOEFF, vecTEDVE, valDENSITY, valNELE, matSPANDIR, vecTE, matPLEX, matWCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG, ...
                matVUINF, matWVLST, vecWLE, vecWLEDVE, matWELST, valAREA, valSPAN, vecDVESYM, vecWDVESYM, vecDGAMMA_DT, vecDGAMMA_DETA, valWSIZE, matWDVE, vecWDVEFLIP, matWEIDX, vecWEMU, vecWVMU, matINTCIRC, matLIFTFREE, matSIDEFREE, matLIFTIND, matSIDEIND, matDRAGIND, matDVEDRAG_DIR, matDVELIFT_DIR, matDVESIDE_DIR, matDVEGRID, vecDVEAREA);
            
            if strcmpi(strATYPE{1}, 'WING')
                [CL(valTIMESTEP), CDi(valTIMESTEP), CY(valTIMESTEP), e(valTIMESTEP)] = fcnWFORCES(valTIMESTEP, vecDVELIFT, vecDVEDRAG, vecDVESIDE, valDENSITY, valAREA, valSPAN);
            else
                [CT(valTIMESTEP), vecDVETHRUST] = fcnPFORCES(strATYPE{1}, valTIMESTEP, vecDVELIFT, vecDVEDRAG, vecDVESIDE, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR, valDENSITY, valDIAM, valRPM);
            end
        end
    end
    
end

% len = size(matLIFTFREE,1);
% for i = 1:len
%     if i <= 2
%         matDGAMMADT(i,:) = (-matINTCIRC(i+2,:) + 4.*matINTCIRC(i+1,:) - 3.*matINTCIRC(i,:))./(2.*valDELTIME);
%     elseif i >= 3 && i <= len-2
%         matDGAMMADT(i,:) = ((-matINTCIRC(i+2,:) + 8.*matINTCIRC(i+1,:) - 8.*matINTCIRC(i-1,:) + matINTCIRC(i-2,:))./(12.*valDELTIME));
%     elseif i == len - 1
%         matDGAMMADT(i,:) = ((matINTCIRC(i+1,:) - matINTCIRC(i,:))./valDELTIME);
%     elseif i == len
%         matDGAMMADT(i,:) = ((matINTCIRC(i,:) - matINTCIRC(i-1,:))./valDELTIME);
%     end
% end

len = size(matLIFTFREE,1);
lambda = 0.5;
matDGAMMADT(1,:) = matINTCIRC(1,:).*0;
for i = 3:len
    if i <= len-2
        matDGAMMADT(i,:) = ((-matINTCIRC(i+2,:) + 8.*matINTCIRC(i+1,:) - 8.*matINTCIRC(i-1,:) + matINTCIRC(i-2,:))./(12.*valDELTIME));
    elseif i == len - 1
        matDGAMMADT(i,:) = ((matINTCIRC(i+1,:) - matINTCIRC(i-1,:))./(2.*valDELTIME));
    elseif i == len
        matDGAMMADT(i,:) = ((matINTCIRC(i,:) - matINTCIRC(i-1,:))./valDELTIME);
    end
end

% len = size(matLIFTFREE,1);
% lambda = 0.5;
% matDGAMMADT(1,:) = matINTCIRC(1,:).*0;
% for i = 2:len
%     matDGAMMADT(i,:) = ((matINTCIRC(i,:) - matINTCIRC(i-1,:))./valDELTIME).*lambda + (1 - lambda).*matDGAMMADT(i-1,:);
% end

tmpLIFTFREE = matLIFTFREE + matDGAMMADT;
tmpDVELIFT = tmpLIFTFREE + matLIFTIND;
tmpDVEDRAG = matDRAGIND;
tmpDVESIDE = matSIDEFREE + matSIDEIND;

for i = 1:size(matLIFTFREE,1)
    tmpDVETHRUST(i,:) = dot(matDVELIFT_DIR(:,:,i).*tmpDVELIFT(i,:)', repmat([0 0 1], length(tmpDVELIFT(i,:)), 1), 2) ...
        + dot(matDVEDRAG_DIR(:,:,i).*tmpDVEDRAG(i,:)', repmat([0 0 1], length(tmpDVEDRAG(i,:)), 1), 2) ...
        + dot(matDVESIDE_DIR(:,:,i).*tmpDVESIDE(i,:)', repmat([0 0 1], length(tmpDVESIDE(i,:)), 1), 2);
    
    if strcmpi(strATYPE, 'PROPELLER')
        CT_U(i,1) = nansum(tmpDVETHRUST(i,:))./(valDENSITY.*((valRPM/60).^2).*(valDIAM.^4));
    else
        CT_U(i,1) = nansum(tmpDVETHRUST(i,:))./(valDENSITY.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
    end
    
end

save('coarse_fixed_2.mat');
figure(20);
plot(CT, '-k')
hold on
plot(CT_U, '--b');
hold off
grid minor
box on
axis tight

% CL2 = sum(tmpLIFTFREE,2)./(0.5.*valDENSITY.*valAREA.*(valUINF^2));
% CLdg = sum(matDGAMMADT,2)./(0.5.*valDENSITY.*valAREA.*(valUINF^2));
% CLc = sum(matINTCIRC,2)./(0.5.*valDENSITY.*valAREA.*(valUINF^2))
% valAR = (valSPAN.^2)./valAREA;
% CL2D = CL2.*((valAR + 2)/valAR);
%
% hFig23 = figure(23);
% clf(23);
% load('Stuff/Gust Response/Kussner.mat')
% plot(Kussner(:,1).*0.5 + (valGUSTSTART+1).*valDELTIME, smooth(Kussner(:,2)), '-k', 'LineWidth', 1);
% box on
% axis tight
% grid minor
%
% hold on
% plot(valDELTIME.*[valGUSTSTART+1:valMAXTIME], CL2D(valGUSTSTART+1:end), '-b', 'LineWidth', 1);
% % plot(valDELTIME.*[valGUSTSTART+1:valMAXTIME], CLdg(valGUSTSTART+1:end), '--r', 'LineWidth', 1.5);
% % plot(valDELTIME.*[valGUSTSTART+1:valMAXTIME], CLc(valGUSTSTART+1:end), '-.k', 'LineWidth', 1.5);
% hold off
%
% xlabel('Distance Travelled by Wing (m)');
% ylabel('C_l')
%
% legend('Kussner Function','C_L','dGammadt','Integrated Circulation','Location','SouthEast')


% hFig24 = figure(24);
% clf(24);
% load('Stuff/Gust Response/ZAERO.mat')
% plot(smooth(ZAERO(:,1)), smooth(ZAERO(:,2)), '-k');
% hold on
%
% load('Stuff/Gust Response/UVLM.mat')
% plot(smooth(UVLM(:,1)), smooth(UVLM(:,2)), '-.k');
%
% plot(0.135+((valDELTIME.*([0:valMAXTIME-1]))), CL2, '--b');
% box on
% axis tight
% grid minor
% %profile viewer
%
% GOLAND_X = 0.135+((valDELTIME.*([0:valMAXTIME-1])));
