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
strFILE = 'inputs/simple_wing.dat';
% strFILE = 'inputs/simple_wing2.dat';
% strFILE = 'inputs/simple_wing_test.dat';

[matPOINTS, strATYPE, vecSYM, flagRELAX, valMAXTIME, valDELTIME, valALPHA, valBETA, matTEPOINTS, matLEPOINTS, vecULS] = fcnOPREAD(strFILE);
[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matELOC, matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG, matCONTROL, vecDVEAREA] = fcnTRIANG(matPOINTS);
[vecLE, vecLEDVE, vecTE, vecTEDVE, matSPANDIR] = fcnLETEGEN(strATYPE, valNELE, matVLST, matELST, matDVECT, matEATT, matLEPOINTS, matTEPOINTS);

flagRELAX = 0;
valMAXTIME = 20
valDENSITY = 1.225;

matUINF = repmat(fcnUINFWING(valALPHA, 0), valNELE, 1);

% [hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], [], matROTANG, [3 1 4 4], 'opengl');
% view([33, 28])

%% D-Matrix Creation
matKINCON_P = matCONTROL;
matKINCON_DVE = [1:valNELE]';

% Points where flow tangency is enforced
matD = fcnDWING9(strATYPE, matEATT, matPLEX, valNELE, matELOC, matELST, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecLEDVE, vecTEDVE, vecSYM, matROTANG, matSPANDIR, matKINCON_P, matKINCON_DVE);
valDLEN = length(matD);

%% Preparing to timestep
hFig21 = []; matWADJE = []; matWELST = []; matWDVE = []; valWNELE = [];
matWEATT = []; matWEIDX = []; matWELOC = []; matWPLEX = []; matWDVECT = [];
matWALIGN = []; matWVATT = []; matWVNORM = []; matWCENTER = [];
matWAKEGEOM = []; matWCOEFF = []; matWVLST = []; matWROTANG = [];

valWSIZE = length(vecTE);

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
    %   Solve entire WD-Matrix (?)
    %   Relaxation procedure (Relax, create W-Matrix and W-Resultant, solve W-Matrix)
    %   Calculate surface normal forces
    %   Calculate DVE normal forces
    %   Calculate induced drag
 
    % Moving the wing
    [matVLST, matCENTER, matNEWWAKE, matCONTROL, matKINCON_P] = fcnMOVEWING(matUINF, valDELTIME, matVLST, matCENTER, matELST, vecTE, matCONTROL, matKINCON_P);
    
    % Generating new wake elements
    if any(vecTE)
        [matWAKEGEOM, matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC,...
            matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWCOEFF, matWROTANG, vecWLE, vecWTE, vecWLEDVE, vecWTEDVE] = fcnCREATEWAKE(valTIMESTEP, strATYPE, vecULS, matNEWWAKE, matWAKEGEOM, matCOEFF, valWSIZE, ...
            vecTE, vecTEDVE, matCENTER, matROTANG, matWCOEFF);
        
        % Rebuild wing resultant
        vecR = fcnRWING(valDLEN, valTIMESTEP, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matKINCON_P, matKINCON_DVE, matDVECT);
        matCOEFF = fcnSOLVED(matD, vecR, valNELE);
        matCOEFF_HSTRY(:,:,valTIMESTEP + 1) = matCOEFF;
        
        % Update wake coefficients
        matWCOEFF(end - valWSIZE*2 + 1:end, :) = fcnDWAKENEW(valTIMESTEP, strATYPE, vecULS, valWNELE, vecWLE, vecWLEDVE, vecWTE, vecWTEDVE, matWEATT, matWELST, matWROTANG, matWCENTER, matWVLST, vecTE, vecTEDVE, matCOEFF, matCENTER, matROTANG, matWCOEFF); 
    end
    
    if flagRELAX == 1 && valTIMESTEP > 10
        [hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], matUINF, matROTANG, [3 1 4 4], 'opengl');    
        [hFig1] = fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER);       
        
        [matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG, matWAKEGEOM] = ...
            fcnRELAX(matUINF, valDELTIME, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, ...
            matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST);
    end

    hFig21 = fcnPLOTCOEFF(hFig21, valTIMESTEP, matCOEFF_HSTRY);
    
    % Forces
    % Jank is in session, The Right Honourable FORLOOP is preciding
    % VECTORIZE IT TRAVIS
    vecDVELIFT = nan(valNELE,1);
    for jj = 1:size(vecTE,1)
        pts = matVLST(matELST(vecTE(jj),:),:); % End points of TE edge segement
        pts_loc = fcnGLOBSTAR(pts - matCENTER(vecTEDVE(jj),:), repmat(matROTANG(vecTEDVE(jj),:),2,1)); % In element local (The TE DVE where this edge is from)
        
        len = 100; % Number of divisions of this line (Jank)
        points = [linspace(pts_loc(1,1), pts_loc(2,1), len)' linspace(pts_loc(1,2), pts_loc(2,2), len)' linspace(pts_loc(1,3), pts_loc(2,3), len)'];
        distance = sqrt(sum((pts_loc(2,:) - pts_loc(1,:)).^2,2)); % Length of entire TE edge
        vec = (pts_loc(2,:) - pts_loc(1,:))./distance; % Direction of this edge (in local)
        
        % Circulation at the points along the edge (oriented along the edge)
        circ = sum([0.5.*points(:,2).^2 points(:,2) 0.5.*points(:,1).^2 points(:,1) ones(size(points(:,1)))].*matCOEFF(vecTEDVE(jj),:),2).*vec; 
        circ = fcnSTARGLOB(circ, repmat(matROTANG(vecTEDVE(jj),:), len, 1)); % Translate to global
        F = cross( repmat(valDENSITY.*matUINF(vecTEDVE(jj),:), len, 1), circ, 2); % A special guest mix on the track, Kutty J
        vecDVELIFT(vecTEDVE(jj),1) = sqrt(sum((sum(F,1).*(distance/len)).^2,2)); % Multiply by the length of the discretization, and sum
    end

    CL = nansum(vecDVELIFT)./(0.5.*valDENSITY.*sum(vecDVEAREA));
    CDi = nan;
    fprintf('Timestep: %d\t\tCL = %0.5f\t\tCDi = %0.5f\n', valTIMESTEP, CL, CDi);
end

% Plot
[hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], matUINF, matROTANG, [3 1 4 4], 'opengl');
[hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), matUINF, matROTANG, 'xr', 20);
if valTIMESTEP > 0
    [hFig1] = fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER);
    [hFig1] = fcnPLOTCIRC(hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, matWPLEX, real(matWCOEFF), matUINF, matWROTANG, 'xb', 4);
end

% s_ind = fcnSDVEVEL(matCENTER, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);
% w_ind = fcnSDVEVEL(matCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER);
% q_ind = s_ind + w_ind + matUINF;
% hold on
% quiver3(matCENTER(:,1), matCENTER(:,2), matCENTER(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3),'m')
% hold off
        
