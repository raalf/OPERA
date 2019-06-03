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
strFILE = 'inputs/ellipse.dat';
% strFILE = 'inputs/test2.dat';

[matPOINTS, strATYPE, vecSYM, flagRELAX, valMAXTIME, valDELTIME, valALPHA, ...
    valBETA, matTEPOINTS, matLEPOINTS, vecULS, valAREA, valSPAN, valDENSITY, vecDVESYM, valDIAM, valCOLL, valRPM, valJ] = fcnOPREAD(strFILE);
[TR, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matELOC, matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG, matCONTROL, vecDVEAREA]...
    = fcnTRIANG(matPOINTS, 'SURFACE', []);
[vecLE, vecLEDVE, vecTE, vecTEDVE, matSPANDIR, vecSYM, vecSYMDVE] = fcnLETEGEN(strATYPE, valNELE, matVLST, matELST, matDVECT, matEATT, matLEPOINTS, matTEPOINTS, vecSYM);

flagRELAX = 1
flagGIF = 1;
valMAXTIME = 60

% matUINF
if strcmpi(strATYPE{1}, 'ROTOR')
    rps = valRPM./60;
    rev = 2*pi*rps*valDELTIME;
    
    dist = 2.*pi.*matCENTER(:,2).*(rev/(2*pi));
    dir = cross(repmat([0 0 1], size(matCENTER,1), 1), matCENTER, 2);
    dir = dir./sqrt(sum(dir.^2, 2));
    matUINF = dir.*(dist./valDELTIME) + valJ.*rps.*valDIAM.*fcnUINFWING(valALPHA, 0);
    
    dist = 2.*pi.*matVLST(:,2).*(rev/(2*pi));
    dir = cross(repmat([0 0 1], size(matVLST,1), 1), matVLST, 2);
    dir = dir./sqrt(sum(dir.^2, 2));
    matVUINF = dir.*(dist./valDELTIME) + valJ.*rps.*valDIAM.*fcnUINFWING(valALPHA, 0);
else
    matUINF = repmat(fcnUINFWING(valALPHA, 0), valNELE, 1);
    matVUINF = repmat(fcnUINFWING(valALPHA, 0), size(matVLST,1), 1);
end

if flagGIF == 1
    if exist('GIF','dir') ~= 7
        mkdir('GIF');
    end
    valGIFNUM = size(dir('GIF'),1) - 1;
end

hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], matUINF, matROTANG, [], 'opengl');

%% D-Matrix Creation
matKINCON_P = matCONTROL;
matKINCON_DVE = [1:valNELE]';

% Points where flow tangency is enforcedfe
matD = fcnDWING9(strATYPE, matEATT, matPLEX, valNELE, matELOC, matELST, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecLEDVE, vecTEDVE, matROTANG, matSPANDIR, matKINCON_P, matKINCON_DVE, vecSYM, vecSYMDVE, vecDVESYM);
valDLEN = length(matD);

%% Preparing to timestep
% Initializing
matWELST = []; matWDVE = []; valWNELE = [];
matWEATT = []; matWEIDX = []; matWELOC = []; matWPLEX = []; matWDVECT = [];
matWALIGN = []; matWVATT = []; matWVNORM = []; matWCENTER = [];
matWAKEGEOM = []; matWCOEFF = []; matWVLST = []; matWROTANG = [];
matWPLANE = [];
vecWDVECIRC = []; CL = nan(valMAXTIME,1); CDi = nan(valMAXTIME,1);
e = nan(valMAXTIME,1); gust_vel_old = matCENTER.*0; vecWDVESYM = logical([]); vecWSYMDVE = []; vecWSYM = [];

valWSIZE = length(vecTE);

% Building wing resultant
vecR = fcnRWING(valDLEN, 0, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matKINCON_P, matKINCON_DVE, matDVECT, []);

% Solving for wing coefficients
matCOEFF = fcnSOLVED(matD, vecR, valNELE);
matCOEFF = fcnADJCOEFF(matVLST, matVATT, matCENTER, matROTANG, matDVE, matCOEFF, matELST, matEATT, matEIDX, vecTE, valNELE, []);
matCOEFF_HSTRY(:,:,1) = matCOEFF;

valGUSTAMP = 0.0025;
valGUSTL = 7.3152;
flagGUSTMODE = 2;
valGUSTSTART = 12;

if strcmpi(strATYPE{1}, 'WING')
    valPRESTEPS = 0
    valDELTIME = valDELTIME*10;
    strWAKE_TYPE = 'STEADY';
else
    valPRESTEPS = 0;
    strWAKE_TYPE = strATYPE{3};
end

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
    
    if strcmpi(strATYPE{1}, 'WING') && valTIMESTEP == valPRESTEPS + 1
        valDELTIME = valDELTIME/10;
        strWAKE_TYPE = strATYPE{3};
    end
    
    if strcmpi(strATYPE{1}, 'ROTOR')
        [matVLST, matCENTER, matNEWWAKE, matCONTROL, matKINCON_P, matUINF, matVUINF] = fcnMOVEROTOR(matUINF, matVUINF, valRPM, valJ, valDIAM, valALPHA, valDELTIME, matVLST, matCENTER, matELST, vecTE, matCONTROL, matKINCON_P);
    else
        % Moving the wing
        [matVLST, matCENTER, matNEWWAKE, matCONTROL, matKINCON_P] = fcnMOVEWING(matUINF, valDELTIME, matVLST, matCENTER, matELST, vecTE, matCONTROL, matKINCON_P);
    end
    
    % Generating new wake elements
    if any(vecTE)
        [matWAKEGEOM, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC,...
            matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWCOEFF, matWROTANG, vecWLE, vecWTE, vecWLEDVE, vecWTEDVE, vecWDVECIRC, vecWSYM, vecWSYMDVE, vecWDVESYM, matWVGRID, vecWOTE] = fcnCREATEWAKE(valTIMESTEP, strATYPE, vecULS, matNEWWAKE, matWAKEGEOM, matCOEFF, valWSIZE, ...
            vecTE, vecTEDVE, matCENTER, matROTANG, matWCOEFF, matWPLEX, vecWDVECIRC, vecWSYMDVE, vecSYMDVE, vecWDVESYM, vecDVESYM, vecWSYM);
        
        % Rebuild wing resultant
        vecR = fcnRWING(valDLEN, valTIMESTEP, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matKINCON_P, matKINCON_DVE, matDVECT, vecWDVESYM);
        matCOEFF = fcnSOLVED(matD, vecR, valNELE);
        matCOEFF = fcnADJCOEFF(matVLST, matVATT, matCENTER, matROTANG, matDVE, matCOEFF, matELST, matEATT, matEIDX, [vecTE; vecSYM], valNELE, []);
        
        % Update wake coefficients
        if strcmpi(strWAKE_TYPE,'STEADY')
            matWCOEFF = fcnDWAKE('STEADY', valTIMESTEP, strATYPE, vecULS, valWNELE, vecWLE, vecWLEDVE, vecWTE, vecWTEDVE, matWEATT, matWELST, matWROTANG, matWCENTER, matWVLST, vecTE, vecTEDVE, matCOEFF, matCENTER, matROTANG, matWCOEFF, matWPLEX, vecWDVECIRC, vecWSYM, vecWSYMDVE, vecWDVESYM, matWVATT, matWDVE);
        else
            [matWCOEFF(end - valWSIZE*2 + 1:end, :), vecWDVECIRC(end - valWSIZE*2 + 1:end, :)] = fcnDWAKENEW(valTIMESTEP, strATYPE, vecULS, valWNELE, vecWLE, vecWLEDVE, vecWTE, vecWTEDVE, matWEATT, matWELST, matWROTANG, matWCENTER, matWVLST, vecTE, vecTEDVE, matCOEFF, matCENTER, matROTANG, matWCOEFF, matWPLEX, vecWSYM, vecWSYMDVE, vecWDVESYM, matWVATT, matWDVE);
        end
        matWCOEFF = fcnADJCOEFF(matWVLST, matWVATT, matWCENTER, matWROTANG, matWDVE, matWCOEFF, matWELST, matWEATT, matWEIDX, [vecWLE; vecWOTE; vecWSYM], valWNELE, []);
        
        % Relaxing Wake
        if flagRELAX == 1 && valTIMESTEP > valPRESTEPS
            %             hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, [], 'opengl');
            %             fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, valWSIZE, valPRESTEPS);
            [matWELST, matWVLST, matWDVE, valWNELE, matWEIDX, matWPLEX, matWDVECT, matWVATT, matWCENTER, matWROTANG, matWAKEGEOM, matWVGRID] = ...
                fcnRELAX5(valTIMESTEP, valDELTIME, valNELE, matCOEFF, matPLEX, valWNELE, matWCOEFF, matWDVE, matWVLST, matWPLEX, valWSIZE, ...
                matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, matWELST, matWVATT, matWEIDX, vecDVESYM, vecWDVESYM, vecWSYM, matWVGRID, valPRESTEPS);
            
            matWCOEFF = fcnDWAKE('STEADY', valTIMESTEP, strATYPE, vecULS, valWNELE, vecWLE, vecWLEDVE, vecWTE, vecWTEDVE, matWEATT, matWELST, matWROTANG, matWCENTER, matWVLST, vecTE, vecTEDVE, matCOEFF, matCENTER, matROTANG, matWCOEFF, matWPLEX, vecWDVECIRC, vecWSYM, vecWSYMDVE, vecWDVESYM, matWVATT, matWDVE);
            matWCOEFF = fcnADJCOEFF(matWVLST, matWVATT, matWCENTER, matWROTANG, matWDVE, matWCOEFF, matWELST, matWEATT, matWEIDX, [vecWLE; vecWOTE; vecWSYM], valWNELE, []);
        end
                
        % Updating coefficient convergence history
        matCOEFF_HSTRY(:,:,valTIMESTEP + 1) = matCOEFF;
        
        if flagGIF == 1 && valTIMESTEP > valPRESTEPS
            hFig1 = fcnGIF(valTIMESTEP - valPRESTEPS, valNELE, matDVE, matVLST, matCENTER, matELST, matDVECT, matPLEX, matCOEFF, matUINF, matROTANG, ...
                valWNELE, matWDVE, matWVLST, matWCENTER, matWELST, matWDVECT, valWSIZE, valPRESTEPS, matWVGRID, valGIFNUM);
        end
    end
    
    %% Plot wake circ
    
    %         hFig20 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, [], 'opengl');
    %         fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, 'r', 100);
    %         tmatWVLST = repmat(matWVLST(1,:), size(matWVLST,1), 1);
    %         for i = 1:size(matWVGRID,1)
    %             tmp1 = matWVLST(matWVGRID(i,:),:);
    %             tmp2 = [tmp1(1,:); tmp1(1:end-1,:)];
    %             tmatWVLST(matWVGRID(i,:), 2) = cumsum(sqrt(sum((tmp2 - tmp1).^2, 2))) + tmatWVLST(matWVGRID(i,:),2);
    %         end
    %
    %         for i = 1:size(matWVGRID,2)
    %             tmp1 = matWVLST(matWVGRID(:,i),:);
    %             tmp2 = [tmp1(1,:); tmp1(1:end-1,:)];
    %             tmatWVLST(matWVGRID(:,i), 1) = cumsum(sqrt(sum((tmp2 - tmp1).^2, 2))) + tmatWVLST(matWVGRID(:,i),1);
    %         end
    %         tmatWAKEGEOM = [];
    %         tmatWAKEGEOM(:,:,1) = tmatWVLST(matWDVE(:,1,1),:);
    %         tmatWAKEGEOM(:,:,2) = tmatWVLST(matWDVE(:,2,1),:);
    %         tmatWAKEGEOM(:,:,3) = tmatWVLST(matWDVE(:,3,1),:);
    %         tmatWETA = nan(valWNELE, 1);
    %         tmatWETA(reshape([1:valWSIZE.*2:valWNELE]' + [1:valWSIZE]-1, [], 1),1) = 1;
    %         tmatWETA(reshape([1:valWSIZE.*2:valWNELE]' + [1:valWSIZE] + valWSIZE - 1, [], 1),1) = 2;
    %         [~, ~, ~, ~, ~, ~, ~, ~, tmatWPLEX, tmatWDVECT, ~, ~, tmatWCENTER, tmatWROTANG, ~] = fcnTRIANG(tmatWAKEGEOM, 'WAKE', tmatWETA);
    %         tmatWCOEFF = fcnDWAKE('STEADY', valTIMESTEP, strATYPE, vecULS, valWNELE, vecWLE, vecWLEDVE, vecWTE, vecWTEDVE, matWEATT, matWELST, tmatWROTANG, tmatWCENTER, tmatWVLST, vecTE, vecTEDVE, matCOEFF, matCENTER, matROTANG, matWCOEFF, matWPLEX, vecWDVECIRC, vecWSYM, vecWSYMDVE, vecWDVESYM, matWVATT, matWDVE);
    %         tmatWCOEFF = fcnADJCOEFF(tmatWVLST, matWVATT, tmatWCENTER, tmatWROTANG, matWDVE, matWCOEFF, matWELST, matWEATT, matWEIDX, [vecWLE; vecWOTE; vecWSYM], valWNELE, []);
    %         fcnPLOTWAKE(1, hFig20, matWDVE, valWNELE, tmatWVLST, matWELST, tmatWDVECT, tmatWCENTER, valWSIZE, 0);
    %         fcnPLOTCIRC(hFig20, matWDVE, valWNELE, tmatWVLST, matWELST, tmatWDVECT, tmatWCENTER, tmatWPLEX, tmatWCOEFF, matUINF, tmatWROTANG, 'r', 200);
    %         grid minor
    %         box on
    %         view([-29 7]);
    
    %% Calculating Forces
    [CL(valTIMESTEP,1), CDi(valTIMESTEP,1), CY(valTIMESTEP,1), e(valTIMESTEP,1), vecDVELIFT, vecDVEDRAG, matDVEDRAG_DIR, matDVELIFT_DIR, matSIDE_DIR] = fcnFORCES(valTIMESTEP, matVLST, matCENTER, matELST, matROTANG, ...
        matUINF, matCOEFF, vecTEDVE, valDENSITY, valNELE, matSPANDIR, vecTE, vecDVEAREA, matPLEX, matWCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matVUINF, matWVLST, vecWLE, vecWLEDVE, matWELST, valAREA, valSPAN, matWDVECT, matDVECT, vecDVESYM, vecWDVESYM);
end

%% Plotting
% try
%     clf(1);
% end
% hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, [], 'opengl');
% fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, 'r', 100);
% if valMAXTIME > 0 && any(vecTE)
%     fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, valWSIZE, valPRESTEPS);
% %     fcnPLOTCIRC(hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, matWPLEX, matWCOEFF, matUINF, matWROTANG, 'r', 5);
% end
%
% q_ind = fcnSDVEVEL(matCENTER, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM) + matUINF;
% hold on
% quiver3(matCENTER(:,1), matCENTER(:,2), matCENTER(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3))
% hold off
% view([153 15])
%
% granularity = 0.005;
% % x = [-0.1:granularity:0.8];
% x = 0.175;
% % x = 0.9;
% y = [-3:granularity:3];
% % y = 0.3
% z = [0.0000];
% [X,Y,Z] = meshgrid(x,y,z);
% fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');
%
% q_ind = fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER, vecDVESYM);
%
% hold on
% quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3))
% hold off
% view([153 15])
%
% figure(2);
% clf(2);
% plot(y, q_ind(:,3), '-ok');
% box on
% grid minor
% axis tight
