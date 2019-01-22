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
% strFILE = 'inputs/simple_wing.dat';
% strFILE = 'inputs/simple_wing2.dat';
% strFILE = 'inputs/simple_wing_test.dat';
strFILE = 'inputs/ellipse.dat';
% strFILE = 'inputs/simple_wing_single.dat';

[matPOINTS, strATYPE, vecSYM, flagRELAX, valMAXTIME, valDELTIME, valALPHA, valBETA, matTEPOINTS, matLEPOINTS, vecULS] = fcnOPREAD(strFILE);
[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matELOC, matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG, matCONTROL, vecDVEAREA] = fcnTRIANG(matPOINTS);
[vecLE, vecLEDVE, vecTE, vecTEDVE, matSPANDIR] = fcnLETEGEN(strATYPE, valNELE, matVLST, matELST, matDVECT, matEATT, matLEPOINTS, matTEPOINTS);

valALPHA = 4
flagRELAX =  0
valMAXTIME = 100
valDENSITY = 1.225;

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
matWADJE = []; matWELST = []; matWDVE = []; valWNELE = [];
matWEATT = []; matWEIDX = []; matWELOC = []; matWPLEX = []; matWDVECT = [];
matWALIGN = []; matWVATT = []; matWVNORM = []; matWCENTER = [];
matWAKEGEOM = []; matWCOEFF = []; matWVLST = []; matWROTANG = [];
vecWDVECIRC = [];

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
            matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWCOEFF, matWROTANG, vecWLE, vecWTE, vecWLEDVE, vecWTEDVE, vecWDVECIRC] = fcnCREATEWAKE(valTIMESTEP, strATYPE, vecULS, matNEWWAKE, matWAKEGEOM, matCOEFF, valWSIZE, ...
            vecTE, vecTEDVE, matCENTER, matROTANG, matWCOEFF, matWPLEX, vecWDVECIRC);
        
        % Rebuild wing resultant
        vecR = fcnRWING(valDLEN, valTIMESTEP, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER, matKINCON_P, matKINCON_DVE, matDVECT);
        matCOEFF = fcnSOLVED(matD, vecR, valNELE);
        
        % Update wake coefficients
        [matWCOEFF(end - valWSIZE*2 + 1:end, :), vecWDVECIRC(end - valWSIZE*2 + 1:end, :)] = fcnDWAKENEW(valTIMESTEP, strATYPE, vecULS, valWNELE, vecWLE, vecWLEDVE, vecWTE, vecWTEDVE, matWEATT, matWELST, matWROTANG, matWCENTER, matWVLST, vecTE, vecTEDVE, matCOEFF, matCENTER, matROTANG, matWCOEFF, matWPLEX);
        
        if flagRELAX == 1 && valTIMESTEP > 1
            [hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], matUINF, matROTANG, [3 1 4 4], 'opengl');
            [hFig1] = fcnPLOTWAKE(1, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER);
            view([-15 12])
            setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera') ;
            
            [matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG, matWAKEGEOM] = ...
                fcnRELAX2(valTIMESTEP, matWELOC, matWEATT, matWEIDX, matUINF, valDELTIME, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, ...
                matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST);
            
            %             [matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG, matWAKEGEOM] = ...
            %                 fcnRELAX(matUINF, valDELTIME, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, ...
            %                 matROTANG, matWROTANG, matCENTER, matWCENTER, vecWLE, vecWTE, matWELST);
            
            [hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], matUINF, matROTANG, [3 1 4 4], 'opengl');
            [hFig1] = fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER);
            view([-15 12])
            
            matWCOEFF = fcnDWAKE(valTIMESTEP, strATYPE, vecULS, valWNELE, vecWLE, vecWLEDVE, vecWTE, vecWTEDVE, matWEATT, matWELST, matWROTANG, matWCENTER, matWVLST, vecTE, vecTEDVE, matCOEFF, matCENTER, matROTANG, matWCOEFF, matWPLEX, vecWDVECIRC);
        end
        matCOEFF_HSTRY(:,:,valTIMESTEP + 1) = matCOEFF;
    end
    
    %     hFig21 = fcnPLOTCOEFF(valTIMESTEP, matCOEFF_HSTRY);
    
    % Forces
    % Jank is in session, The Right Honourable FORLOOP is preciding
    % VECTORIZE IT TRAVIS
    vecDVELIFT = nan(valNELE,1);
    vecDVEDRAG = nan(valNELE,1);
    matDVEDRAG_DIR = matUINF./sqrt(sum(matUINF.^2,2));
    matDVELIFT_DIR = cross(matDVEDRAG_DIR, matSPANDIR, 2);
    matDVESIDEF_DIR = cross(matDVELIFT_DIR, matDVEDRAG_DIR, 2);
    
    liftfree = nan(valNELE,1);
    liftind = nan(valNELE,1);
    sidefree = nan(valNELE,1);
    sideind = nan(valNELE,1);
    dragind = nan(valNELE,1);
    
    [hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], matUINF, matROTANG, [3 1 4 4], 'opengl');
    [hFig1] = fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER);
    
    tmp1 = [];
    tmp2 = [];
    tmp3 = [];
    for jj = 1:size(vecTE,1)
        % Lift
        pts = matWVLST(matWELST(vecWLE(jj),:),:); % End points of TE edge segement
        pts_loc = fcnGLOBSTAR(pts - matWCENTER(vecWLEDVE(jj),:), repmat(matWROTANG(vecWLEDVE(jj),:),2,1)); % In element local (The TE DVE where this edge is from)
        
        len = 10; % Number of divisions of this line (Jank)
        points = [linspace(pts_loc(1,1), pts_loc(2,1), len)' linspace(pts_loc(1,2), pts_loc(2,2), len)' linspace(pts_loc(1,3), pts_loc(2,3), len)'];
        distance = sqrt(sum((pts_loc(2,:) - pts_loc(1,:)).^2,2)); % Length of entire TE edge
        vec = (pts_loc(1,:) - pts_loc(2,:))./distance; % Direction of this edge (in local)
        
        % Circulation at the points along the edge (oriented along the edge)
        circ = sum([0.5.*points(:,2).^2 points(:,2) 0.5.*points(:,1).^2 points(:,1) ones(size(points(:,1)))].*matCOEFF(vecTEDVE(jj),:),2).*vec;
        circ = fcnSTARGLOB(circ, repmat(matWROTANG(vecWLEDVE(jj),:), len, 1)); % Translate to global
        F = cross( repmat(valDENSITY.*matUINF(vecTEDVE(jj),:), len, 1), circ, 2); % A special guest mix on the track, Kutty J
        spans = repmat(distance/len, len, 1);
        spans([1 len],:) = spans([1 len],:).*0.5;
        liftfree(vecTEDVE(jj),1) = sqrt(sum((sum(F.*spans,1)).^2,2)); % Multiply by the length of the discretization, and sum
        
        % Drag
        fpg = fcnSTARGLOB(points, repmat(matWROTANG(vecWLEDVE(jj),:),len,1)) + matWCENTER(vecWLEDVE(jj),:);
%         fpg = fpg - repmat(matWDVECT(vecWLEDVE(jj),:,3),len,1).*1e-1;
%         w_ind = fcnSDVEVEL(fpg, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER) + fcnSDVEVEL(fpg, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);
        w_ind = fcnSDVEVEL(fpg, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER);
        
        w_ind(1,:) = w_ind(1,:).*0;
        w_ind(end,:) = w_ind(end,:).*0;
        
        tmp1 = [tmp1; fpg];
        tmp2 = [tmp2; w_ind];
        
        F_ind = cross( w_ind, circ, 2);
        spans(2,:) = spans(2,:) + 0.5.*(distance/len);
        spans(end-1,:) = spans(end-1,:) + 0.5.*(distance/len);
        force = sum(F_ind.*spans,1); 
        
        tmp3 = [tmp3; F_ind];
        
        dragind(vecTEDVE(jj),1) = dot(force, matDVEDRAG_DIR(vecTEDVE(jj),:), 2);
        liftind(vecTEDVE(jj),1) = dot(force, matDVELIFT_DIR(vecTEDVE(jj),:), 2);
        
        vecDVEDRAG(vecTEDVE(jj),1) = dragind(vecTEDVE(jj),1);
        vecDVELIFT(vecTEDVE(jj),1) = liftfree(vecTEDVE(jj),1) + liftind(vecTEDVE(jj),1);
        %                   vecDVELIFT(vecTEDVE(jj),1) = liftfree(vecTEDVE(jj),1);
    end
    
    hold on
    quiver3(tmp1(:,1), tmp1(:,2), tmp1(:,3), tmp2(:,1), tmp2(:,2), tmp2(:,3),'b')
    quiver3(tmp1(:,1), tmp1(:,2), tmp1(:,3), tmp3(:,1), tmp3(:,2), tmp3(:,3),'r')
    setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera') ;
    hold off
    view([36 15]);
    
    CL = nansum(vecDVELIFT)./(0.5.*valDENSITY.*sum(vecDVEAREA));
    CDi = nansum(vecDVEDRAG)./(0.5.*valDENSITY.*sum(vecDVEAREA));
    e = (CL.^2)./(pi.*((max(matVLST(:,2))-min(matVLST(:,2))).^2)./sum(vecDVEAREA).*CDi);
    fprintf('Timestep: %d\t\tCL = %0.5f\t\tCDi = %0.5f\t\te = %0.5f\n', valTIMESTEP, CL, CDi, e);
end

% Plot
[hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], matUINF, matROTANG, [3 1 4 4], 'opengl');
[hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), matUINF, matROTANG, 'xr', 20);
if valTIMESTEP > 0
    [hFig1] = fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER);
    [hFig1] = fcnPLOTCIRC(hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, matWPLEX, real(matWCOEFF), matUINF, matWROTANG, 'xb', 4);
end
setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera') ;

% s_ind = fcnSDVEVEL(matCENTER, valNELE, matCOEFF, matPLEX, matROTANG, matCENTER);
% w_ind = fcnSDVEVEL(matCENTER, valWNELE, matWCOEFF, matWPLEX, matWROTANG, matWCENTER);
% q_ind = s_ind + w_ind + matUINF;
% hold on
% quiver3(matCENTER(:,1), matCENTER(:,2), matCENTER(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3),'m')
% hold off

