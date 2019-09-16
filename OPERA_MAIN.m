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
strFILE = 'inputs/ellipse.dat';
% strFILE = 'inputs/goland_wing.dat'
% strFILE = 'inputs/kussner.dat'
% strFILE = 'inputs/box_wing.dat'
% strFILE = 'inputs/TMotor.dat'
% strFILE = 'inputs/TMotor_coarse.dat'
% strFILE = 'inputs/TMotor_nocamber.dat'
% strFILE = 'inputs/Leishman_Rotor.dat'
% strFILE = 'inputs/Caradonna_Rotor.dat'
% strFILE = 'Stuff/Columbia/Columbia_234_Rotor.dat'
% strFILE = 'inputs/test.dat'

[matPOINTS, strATYPE, vecSYM, flagRELAX, valMAXTIME, valDELTIME, valALPHA, ...
    valBETA, matTEPOINTS, matLEPOINTS, vecULS, valAREA, valSPAN, valDENSITY, vecDVESYM, valDIAM, ...
    valCOLL, valRPM, valJ, vecDVESURFACE, vecDVEFLIP, valBLADES, valGUSTAMP, valGUSTL, valGUSTMODE, valGUSTSTART, ...
    strHTYPE, matHINGELOC, matHINGEDIR, matHINGEANG] = fcnOPREAD(strFILE);

[TR, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matPLEX, matDVECT, matVATT, ~, matCENTER, matROTANG, ~, vecDVEAREA, matSPANDIR]...
    = fcnTRIANG(matPOINTS, vecDVEFLIP);

flagRELAX = 0
flagGIF = 1;
flagHVRMOD = true;

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
% matHINGEANG = repmat([0 -6.5 0], 3, 1);


%% Preliminary Geometry Stuff
% Duplicating blades for rotors
if valBLADES > 1
    [valNELE, matVLST, matELST, matDVE, matCENTER, matPLEX, vecDVEAREA, ...
        matEATT, matEIDX, vecDVESURFACE, vecDVEFLIP, vecDVESYM, matDVECT, ...
        matROTANG, matVATT, matTEPOINTS, matLEPOINTS, matSPANDIR, matHINGELOC, matHINGEDIR] = fcnDUPBLADE(valBLADES, valNELE, matVLST, matELST, matDVE, ...
        matCENTER, matPLEX, vecDVEAREA, matEATT, matEIDX, vecDVESURFACE, vecDVEFLIP, vecDVESYM, matDVECT, matROTANG, matVATT, matTEPOINTS, matLEPOINTS, matSPANDIR, matHINGELOC, matHINGEDIR);
end

% Generating leading and trailing edge information
[vecLE, vecLEDVE, vecTE, vecTEDVE, vecSYM, vecSYMDVE, matELST] = fcnLETEGEN(strATYPE, matVLST, matELST, matEATT, matLEPOINTS, matTEPOINTS, vecSYM);

if ~isempty(strHTYPE)
    % Setting blade angles for propellers
    [matHINGELOC, matHINGEDIR, matVLST, matCENTER, matPLEX, matDVECT, matROTANG, matSPANDIR] = fcnSETHANG(matHINGEANG, matHINGELOC, matHINGEDIR, matVLST, matDVE, vecDVESURFACE, valBLADES, vecDVEFLIP, matSPANDIR);
end

% Calculating matUINF for rotor or wing
if strcmpi(strATYPE{1}, 'ROTOR')
    vecHUB = [0 0 0];
    vecROTORRADPS = valRPM.*2.*pi./60;
    translation = valJ.*(valRPM.*(pi/30)).*(valDIAM/2).*fcnUINFWING(valALPHA, 0);
    matUINF = cross(repmat([0,0,-vecROTORRADPS],length(matCENTER(:,1)),1),matCENTER) - translation;
    matVUINF = cross(repmat([0,0,-vecROTORRADPS],length(matVLST(:,1)),1),matVLST) - translation;
elseif strcmpi(strATYPE{1}, 'PROPELLER')
    vecHUB = [0 0 0];
    vecROTORRADPS = valRPM.*2.*pi./60;
    translation = valJ.*(valRPM/60).*(valDIAM).*fcnUINFWING(valALPHA, 0);
    matUINF = cross(repmat([0,0,-vecROTORRADPS],length(matCENTER(:,1)),1),matCENTER) - translation;
    matVUINF = cross(repmat([0,0,-vecROTORRADPS],length(matVLST(:,1)),1),matVLST) - translation;
else
    matUINF = repmat(fcnUINFWING(valALPHA, 0), valNELE, 1);
    matVUINF = repmat(fcnUINFWING(valALPHA, 0), size(matVLST,1), 1);
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


hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], matUINF, matROTANG, [], 'opengl');
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
vecDGAMMA_DT = zeros(size(vecDVESURFACE)); valPRESTEPS = 0; strWAKE_TYPE = strATYPE{3};

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
    
    if ~isempty(valGUSTSTART) && valGUSTSTART > 0 && valTIMESTEP > valPRESTEPS
        tmatUINF = matUINF;
        [matUINF, gust_vel_old] = fcnGUSTWING(matUINF, valGUSTAMP, valGUSTL, valGUSTMODE, valDELTIME, 1, valGUSTSTART, matCENTER, gust_vel_old, valPRESTEPS*valDELTIME*10);
        matAINF = (matUINF - tmatUINF)./valDELTIME;
    end
    
    if (strcmpi(strATYPE{1}, 'ROTOR') || strcmpi(strATYPE{1}, 'PROPELLER'))
        [matVLST, matCENTER, matNEWWAKE, matUINF, matVUINF, matPLEX, matDVECT, matROTANG, vecHUB, matSPANDIR] = fcnMOVEROTOR(valRPM, valJ, valDIAM, valALPHA, valDELTIME, matVLST, matELST, vecTE, matDVE, matCENTER, vecDVEFLIP, vecHUB, matSPANDIR);
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
                %                 matWCOEFF = fcnADJCOEFF(vecWVMU, vecWEMU, matWVLST, matWCENTER, matWROTANG, matWDVE, matWCOEFF, matWELST, matWEIDX, valWNELE);
                
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
        end
        
        %% Calculating Forces
        %     hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, [], 'opengl');
        %     hFig1 = fcnPLOTWAKE(1, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, valWSIZE, valPRESTEPS, matWVGRID, vecWDVESURFACE);
        
        
        
        
        
        
        
        
        
        
        
        
        
%         locations = linspace(0.2, 0.8, 3)';
%         
%         wind = nan(size(locations,1), 3, size(vecWLEDVE,1));
%         fpg_og = nan(size(locations,1), 3, size(vecWLEDVE,1));
%         
%         boundind = false(valWNELE,1);
%         boundind(vecWLEDVE) = true;
%         
%         % Need to adjust wake LE elements to straight, one by one
%         for i = 1:valWSIZE
%             tmpWVLST = matWVLST;
%             
%             tmp_verts(1) = matWDVE(vecWLEDVE(i),2);
%             tmp_verts(2) = matWDVE(vecWLEDVE(i),3);
%             
%             tmp_mp = mean(matWVLST(tmp_verts,:),1);
%             tmp_spandir = matSPANDIR(vecTEDVE(i),:);
%             
%             vec_one = matWVLST(tmp_verts(1),:) - tmp_mp;
%             vec_two = matWVLST(tmp_verts(2),:) - tmp_mp;
%             
%             vert_one = ((dot(vec_one, tmp_spandir, 2)).*tmp_spandir) + tmp_mp;
%             vert_two = ((dot(vec_two, tmp_spandir, 2)).*tmp_spandir) + tmp_mp;
%             
%             % Making tmp variables
%             tmpWVLST(tmp_verts(1),:) = vert_one;
%             tmpWVLST(tmp_verts(2),:) = vert_two;
%             
%             tmpWCENTER = (tmpWVLST(matWDVE(:,1),:) + tmpWVLST(matWDVE(:,2),:) + tmpWVLST(matWDVE(:,3),:))./3;
%             
%             % matWPLEX, matWDVECT, matWROTANG
%             P = permute(reshape(tmpWVLST(matWDVE(:,:)',:)', 3, 3, []), [2 1 3]);
%             DNORM = cross(tmpWVLST(matWDVE(:,2),:) - tmpWVLST(matWDVE(:,3),:), tmpWVLST(matWDVE(:,1),:) - tmpWVLST(matWDVE(:,3),:), 2);
%             DNORM = DNORM./sqrt(sum(DNORM.^2,2));
%             DNORM(vecWDVEFLIP,:) = DNORM(vecWDVEFLIP,:).*-1;
%             [tmpWPLEX, tmpWDVECT, tmpWROTANG] = fcnTRITOLEX(P, DNORM, tmpWCENTER);
%             tmpWCOEFF = fcnADJCOEFF(vecWVMU, vecWEMU, tmpWVLST, tmpWCENTER, tmpWROTANG, matWDVE, matWCOEFF, matWELST, matWEIDX, valWNELE);
%             
%             fpg_og(:,:,i) = locations.*vert_one + (1 - locations).*vert_two;
%             wind(:,:,i) = fcnSDVEVEL(fpg_og(:,:,i), valWNELE, tmpWCOEFF, tmpWPLEX, tmpWROTANG, tmpWCENTER, vecWDVESYM, boundind, 1e-5);
%             
%             clf(1)
%             fcnPLOTWAKE(0, gcf, matWDVE, valWNELE, tmpWVLST, matWELST, tmpWDVECT, tmpWCENTER, valWSIZE, valPRESTEPS, matWVGRID, vecWDVESURFACE);
%             axis equal
%             axis tight
%             grid minor
%             box on
%             
%             hold on
%             quiver3(fpg_og(:,1,i), fpg_og(:,2,i), fpg_og(:,3,i), wind(:,1,i), wind(:,2,i), wind(:,3,i))
%             hold off
%             view([-50 17])
%             
%         end
%         
% 
%         hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, [], 'opengl');
%         hFig1 = fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, valWSIZE, valPRESTEPS, matWVGRID, vecWDVESURFACE);
%         fpg2 = reshape(permute(fpg_og, [2 1 3]), size(fpg_og, 2), [])';
%         wind2 = reshape(permute(wind, [2 1 3]), size(wind, 2), [])';
%         hold on
%         quiver3(fpg2(:,1), fpg2(:,2), fpg2(:,3), wind2(:,1), wind2(:,2), wind2(:,3))
%         hold off
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        [vecDVELIFT, vecDVEDRAG, vecDVESIDE, matDVEDRAG_DIR, matDVELIFT_DIR, matDVESIDE_DIR, vecDGAMMA_DT, vecDGAMMA_DETA] = fcnDVEFORCES(valTIMESTEP, strWAKE_TYPE, matVLST, matELST, matROTANG, ...
            matUINF, matCOEFF, vecTEDVE, valDENSITY, valNELE, matSPANDIR, vecTE, matPLEX, matWCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG, ...
            matVUINF, matWVLST, vecWLE, vecWLEDVE, matWELST, valAREA, valSPAN, vecDVESYM, vecWDVESYM, vecDGAMMA_DT, vecDGAMMA_DETA, valWSIZE, matWDVE, vecWDVEFLIP, matWEIDX, vecWEMU, vecWVMU);
        
        if strcmpi(strATYPE{1}, 'WING')
            [CL(valTIMESTEP), CDi(valTIMESTEP), CY(valTIMESTEP), e(valTIMESTEP)] = fcnWFORCES(valTIMESTEP, vecDVELIFT, vecDVEDRAG, vecDVESIDE, valDENSITY, valAREA, valSPAN);
        else
            [CT(valTIMESTEP), vecDVETHRUST] = fcnPFORCES(strATYPE{1}, valTIMESTEP, vecDVELIFT, vecDVEDRAG, vecDVESIDE, matDVELIFT_DIR, matDVEDRAG_DIR, matDVESIDE_DIR, valDENSITY, valDIAM, valRPM);
        end
    end
end
%profile viewer
