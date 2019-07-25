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
% strFILE = 'inputs/test2.dat';
% strFILE = 'inputs/test1.dat';
% strFILE = 'inputs/test.dat';
% strFILE = 'inputs/goland_wing.dat'
% strFILE = 'inputs/rotor_test.dat'
% strFILE = 'inputs/box_wing.dat'
% strFILE = 'inputs/TMotor.dat'
strFILE = 'inputs/TMotor_coarse.dat'

[matPOINTS, strATYPE, vecSYM, flagRELAX, valMAXTIME, valDELTIME, valALPHA, ...
    valBETA, matTEPOINTS, matLEPOINTS, vecULS, valAREA, valSPAN, valDENSITY, vecDVESYM, valDIAM, valCOLL, valRPM, valJ, vecDVESURFACE, vecDVEFLIP, valBLADES] = fcnOPREAD(strFILE);
[TR, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matPLEX, matDVECT, matVATT, ~, matCENTER, matROTANG, ~, vecDVEAREA]...
    = fcnTRIANG(matPOINTS, vecDVEFLIP);

if valBLADES > 1
    [valNELE, matVLST, matELST, matDVE, matCENTER, matPLEX, vecDVEAREA, ...
        matEATT, matEIDX, vecDVESURFACE, vecDVEFLIP, vecDVESYM, matDVECT, ...
        matROTANG, matVATT, matTEPOINTS, matLEPOINTS] = fcnDUPBLADE(valBLADES, valNELE, matVLST, matELST, matDVE, ...
        matCENTER, matPLEX, vecDVEAREA, matEATT, matEIDX, vecDVESURFACE, vecDVEFLIP, vecDVESYM, matDVECT, matROTANG, matVATT, matTEPOINTS, matLEPOINTS);
end

[vecLE, vecLEDVE, vecTE, vecTEDVE, matSPANDIR, vecSYM, vecSYMDVE, matELST] = fcnLETEGEN(strATYPE, valNELE, matVLST, matELST, matDVECT, matEATT, matLEPOINTS, matTEPOINTS, vecSYM);

flagRELAX = 0
flagGIF = 1;

%%
% matUINF
if strcmpi(strATYPE{1}, 'ROTOR')
    vecHUB = [0 0 0];
    vecROTORRADPS = valRPM.*2.*pi./60;
    translation = valJ.*(valRPM./60).*valDIAM.*fcnUINFWING(valALPHA, 0);
    matUINF = cross(repmat([0,0,-vecROTORRADPS],length(matCENTER(:,1)),1),matCENTER) + translation;
    matVUINF = cross(repmat([0,0,-vecROTORRADPS],length(matVLST(:,1)),1),matVLST) + translation;
else
    matUINF = repmat(fcnUINFWING(valALPHA, 0), valNELE, 1);
    matVUINF = repmat(fcnUINFWING(valALPHA, 0), size(matVLST,1), 1);
end

if flagGIF == 1
    if exist('GIF','dir') ~= 7
        mkdir('GIF');
    end
    clear dir
    valGIFNUM = size(dir('GIF'),1) - 1;
end

% hFig1 = fcnPLOTBODY(1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], matUINF, matROTANG, [], 'opengl');

%% D-Matrix Creation
matKINCON_DVE = [1:valNELE]';

% Points where flow tangency is enforcedfe
matD = fcnDWING9(strATYPE, matEATT, matPLEX, valNELE, matELST, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecLEDVE, vecTEDVE, matROTANG, matSPANDIR, matCENTER, matKINCON_DVE, vecSYM, vecSYMDVE, vecDVESYM);
valDLEN = length(matD);

%% Preparing to timestep
% Initializing
matWELST = []; matWDVE = []; valWNELE = [];
matWEATT = []; matWEIDX = []; matWPLEX = []; matWDVECT = [];
matWALIGN = []; matWCENTER = [];
matWAKEGEOM = []; matWCOEFF = []; matWVLST = []; matWROTANG = [];
matWPLANE = []; vecWVMU = []; vecWEMU = []; matWVGRID = [];
matWEGRID = []; vecWDVEFLIP = logical([]); vecWLEDVE = []; vecWTEDVE = [];
vecWDVECIRC = []; CL = nan(valMAXTIME,1); CDi = nan(valMAXTIME,1); CT = nan(valMAXTIME,1);

e = nan(valMAXTIME,1); gust_vel_old = matCENTER.*0; vecWDVESYM = logical([]); vecWSYMDVE = []; vecWSYM = [];

valWSIZE = length(vecTE);

% Building wing resultant
vecR = fcnRWING(valDLEN, 0, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matCENTER, matKINCON_DVE, matDVECT, []);

% Solving for wing coefficients
matCOEFF = fcnSOLVED(matD, vecR, valNELE);
[vecVMU, vecEMU] = fcnVEMU(matVLST, matVATT, matCENTER, matROTANG, matCOEFF, matELST, matEATT, vecTE);
matCOEFF = fcnADJCOEFF(vecVMU, vecEMU, matVLST, matCENTER, matROTANG, matDVE, matCOEFF, matELST, matEIDX, valNELE);
matCOEFF_HSTRY(:,:,1) = matCOEFF;

hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, [], 'opengl');
fcnPLOTCIRC(gcf, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, 'b', 30);
% view([90 90])

% valGUSTAMP = 0.0025;
% valGUSTL = 7.3152;
valGUSTAMP = 0.10;
valGUSTL = 4;
flagGUSTMODE = 2;
valGUSTSTART = 0;

if strcmpi(strATYPE{1}, 'WING')
    valPRESTEPS = 20
    valDELTIME = valDELTIME*10;
    strWAKE_TYPE = 'STEADY';
else
    valPRESTEPS = 0
    strWAKE_TYPE = strATYPE{3};
end

%% Timestep to solution
for valTIMESTEP = 1:valMAXTIME
    if strcmpi(strATYPE{1}, 'WING') && valTIMESTEP == valPRESTEPS + 1
        valDELTIME = valDELTIME/10;
        strWAKE_TYPE = strATYPE{3};
    end
    
    if valGUSTSTART > 0 && valTIMESTEP > valPRESTEPS
        [matUINF, gust_vel_old] = fcnGUSTWING(matUINF, valGUSTAMP, valGUSTL, flagGUSTMODE, valDELTIME, 1, valGUSTSTART, matCENTER, gust_vel_old, valPRESTEPS*valDELTIME*10);
    end
    
    if strcmpi(strATYPE{1}, 'ROTOR')
        [matVLST, matCENTER, matNEWWAKE, matUINF, matVUINF, matPLEX, matDVECT, matROTANG, vecHUB] = fcnMOVEROTOR(matUINF, matVUINF, valRPM, valJ, valDIAM, valALPHA, valDELTIME, matVLST, matELST, vecTE, matDVE, matCENTER, vecDVEFLIP, vecHUB);
    else
        [matVLST, matCENTER, matNEWWAKE] = fcnMOVEWING(matUINF, valDELTIME, matVLST, matCENTER, matELST, vecTE);
    end
    
    % Generating new wake elements
    if any(vecTE)
        [matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWPLEX, matWDVECT, matWCENTER, matWCOEFF, matWROTANG, vecWLE, vecWTE, ...
            vecWLEDVE, vecWTEDVE, vecWDVECIRC, vecWSYM, vecWSYMDVE, vecWDVESYM, matWVGRID, vecWOTE, ...
            vecWOTEDVE, matWEGRID, matWE2GRID, vecWVMU, vecWEMU, vecWDVEFLIP] = fcnCREATEWAKE2(valTIMESTEP, strATYPE, matNEWWAKE, matCOEFF, valWSIZE, ...
            vecTEDVE, matCENTER, matROTANG, matWCOEFF, matWPLEX, vecWDVECIRC, vecWSYMDVE, vecSYMDVE, vecWDVESYM, vecDVESYM, ...
            vecWSYM, matWVGRID, matWVLST, matWELST, matWEGRID, vecWVMU, vecWEMU, vecWDVEFLIP, matWCENTER, matWROTANG, matWDVECT, ...
            matWDVE, matWEIDX, vecWLEDVE, matWEATT);
        
        % Rebuild wing resultant
        vecR = fcnRWING(valDLEN, valTIMESTEP, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matCENTER, matKINCON_DVE, matDVECT, vecWDVESYM);
        matCOEFF = fcnSOLVED(matD, vecR, valNELE);
        [vecVMU, vecEMU] = fcnVEMU(matVLST, matVATT, matCENTER, matROTANG, matCOEFF, matELST, matEATT, vecTE);
        matCOEFF = fcnADJCOEFF(vecVMU, vecEMU, matVLST, matCENTER, matROTANG, matDVE, matCOEFF, matELST, matEIDX, valNELE);
        
        % Update wake coefficients
        [vecWVMU, vecWEMU] = fcnWAKEMU(strATYPE, vecWLE, matWVGRID, matWEGRID, matWE2GRID, vecWVMU, vecWEMU, matWELST, matWVLST, vecTEDVE, matCOEFF, matCENTER, matROTANG);
        matWCOEFF = fcnADJCOEFF(vecWVMU, vecWEMU, matWVLST, matWCENTER, matWROTANG, matWDVE, matWCOEFF, matWELST, matWEIDX, valWNELE);
        
        % Relaxing Wake
        if flagRELAX == 1 && valTIMESTEP > valPRESTEPS
            %             hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, [], 'opengl');
            %             fcnPLOTWAKE(0, gcf, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, valWSIZE, valPRESTEPS, matWVGRID);
            [matWELST, matWVLST, matWDVE, valWNELE, matWEIDX, matWPLEX, matWDVECT, matWCENTER, matWROTANG, matWVGRID] = ...
                fcnRELAX7(valDELTIME, valNELE, matCOEFF, matPLEX, valWNELE, matWCOEFF, matWDVE, matWVLST, matWPLEX, valWSIZE, ...
                matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST, matWEIDX, vecDVESYM, vecWDVESYM, vecWSYM, matWVGRID, vecWDVEFLIP, valPRESTEPS, matWE2GRID);
            
            % Update all wake coefficients
            matWCOEFF = fcnADJCOEFF(vecWVMU, vecWEMU, matWVLST, matWCENTER, matWROTANG, matWDVE, matWCOEFF, matWELST, matWEIDX, valWNELE);
        end
        
        % Updating coefficient convergence history
        matCOEFF_HSTRY(:,:,valTIMESTEP + 1) = matCOEFF;
        
        if flagGIF == 1 && valTIMESTEP > valPRESTEPS
            hFig1 = fcnGIF(valTIMESTEP - valPRESTEPS, valNELE, matDVE, matVLST, matCENTER, matELST, matDVECT, matPLEX, matCOEFF, matUINF, matROTANG, ...
                valWNELE, matWDVE, matWVLST, matWCENTER, matWELST, matWDVECT, valWSIZE, valPRESTEPS, matWVGRID, valGIFNUM);
        end
    end
    
    %     hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, [], 'opengl');
    %     fcnPLOTWAKE(0, gcf, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, valWSIZE, valPRESTEPS, matWVGRID);
    %     fcnPLOTCIRC(gcf, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, matWPLEX, matWCOEFF, matUINF, matWROTANG, 'b', 100);
    %     fcnPLOTCIRC(gcf, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, 'b', 100);
    
    %% Calculating Forces
    if strcmpi(strATYPE{1}, 'WING')
        [CL(valTIMESTEP,1), CDi(valTIMESTEP,1), CY(valTIMESTEP,1), e(valTIMESTEP,1), vecDVELIFT, vecDVEDRAG, matDVEDRAG_DIR, matDVELIFT_DIR, matSIDE_DIR] = fcnWFORCES(valTIMESTEP, matVLST, matCENTER, matELST, matROTANG, ...
            matUINF, matCOEFF, vecTEDVE, valDENSITY, valNELE, matSPANDIR, vecTE, vecDVEAREA, matPLEX, matWCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matVUINF, matWVLST, vecWLE, vecWLEDVE, matWELST, valAREA, valSPAN, matWDVECT, matDVECT, vecDVESYM, vecWDVESYM);
    elseif strcmpi(strATYPE{1}, 'ROTOR')
        CT(valTIMESTEP,1) = fcnRFORCES(valTIMESTEP, matVLST, matCENTER, matELST, matROTANG, ...
            matUINF, matCOEFF, vecTEDVE, valDENSITY, valNELE, matSPANDIR, vecTE, vecDVEAREA, matPLEX, matWCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matVUINF, matWVLST, vecWLE, vecWLEDVE, matWELST, valAREA, valSPAN, matWDVECT, matDVECT, vecDVESYM, vecWDVESYM, valDIAM, valRPM);
    end
    
    %         hFig21 = figure(21);
    %         clf(21);
    %         count = 1;
    %         for jj = 1:2:size(matWVGRID,2).*2
    %             hold on
    %             plot3(1:size(matWVGRID(:,count),1)', repmat(jj, size(matWVGRID,1), 1), vecWVMU(matWVGRID(:,count)), '-ok');
    %             hold off
    %             count = count + 1;
    %         end
    %
    %         count = 1;
    %         for jj = 2:2:size(matWEGRID,2).*2
    %             hold on
    %             plot3(linspace(1,size(matWVGRID(:,count),1),size(matWEGRID(:,count),1)), repmat(jj, size(matWEGRID,1), 1), vecWEMU(matWEGRID(:,count)), '--^m');
    %             hold off
    %             count = count + 1;
    %         end
    %
    %         count = 1;
    %         for jj = 1:2:size(matWE2GRID,2).*2
    %             hold on
    %             plot3(1.5:1:(size(matWE2GRID(:,count),1)+0.5)', repmat(jj, size(matWE2GRID,1), 1), vecWEMU(matWE2GRID(:,count)), '-.ks');
    %             hold off
    %             count = count + 1;
    %         end
    %
    %         grid minor
    %         box on
    %         axis tight
    %         axis equal
    
end

% v1 = matELST(vecLE(1),:);
% v2 = matELST(vecLE(2),:);
%
% len = 200;
% x = repmat(0.5,len,1);
% y = [linspace(matVLST(v1(1),2), matVLST(v1(2),2), len/2)'; linspace(matVLST(v2(1),2), matVLST(v2(2),2), len/2)'];
% z = [linspace(matVLST(v1(1),3), matVLST(v1(2),3), len/2)'; linspace(matVLST(v2(1),3), matVLST(v2(2),3), len/2)'];
% fpg = [x y z];
%
% % fpg = [0.5 0.5 0];
% q_ind = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM,[], 0);
% q_ind2 = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM,[], 1e-20);
%
% hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, [], 'opengl');
% hold on
% quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3))
% hold off
% view([153 15])
%
% p1 = sqrt(sum(q_ind.^2,2));
% p2 = sqrt(sum(q_ind2.^2,2));
% figure(2);
% clf(2);
% hold on
% plot(y, p1, '-ok');
% plot(y, p2, '--^m');
% hold off
% box on
% grid minor
% axis tight
