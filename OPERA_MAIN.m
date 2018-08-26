clear
clc
tic
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
%
strFILE = 'inputs/simple_wing2.dat';
% strFILE = 'inputs/KORDY30.dat';
% strFILE = 'inputs/standard_cirrus.dat';
% strFILE = 'inputs/2dve.dat';
% strFILE = 'inputs/4dve.dat';
% strFILE = 'inputs/4dve_nosym.dat';
% strFILE = 'inputs/nonplanar.dat';

[matPOINTS, strATYPE, vecSYM, flagRELAX, valMAXTIME, valDELTIME, valALPHA, valBETA, matTEPOINTS, matLEPOINTS] = fcnOPREAD(strFILE);

[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matELOC, matPLEX, matDVECT, matVATT, matVNORM, matCENTER, matROTANG] = fcnTRIANG(matPOINTS);

flagRELAX = 0
valMAXTIME = 0
valDENSITY = 1.225

vecUINF = fcnUINFWING(valALPHA, 0);
matUINF = repmat(vecUINF,valNELE,1);

if ~isempty(matTEPOINTS) && ~isempty(matLEPOINTS)
    [vecTE, vecLE] = fcnTELE(matTEPOINTS, matLEPOINTS, matVLST, matELST);
else
    vecTE = [];
    vecLE = [];
end

[hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], matUINF, matROTANG, [3 1 4 4], 'opengl')
view([33, 28])

%% D-Matrix Creation
vecTEDVE = [];
vecSPANDIR = [];

vecLEDVE = nonzeros(sort(matEATT(vecLE,:),2,'descend'));
if ~isempty(vecTE)
    vecTEDVE = nonzeros(sort(matEATT(vecTE,:),2,'descend')); % A vector of trailing edge HDVEs, which corresponds to vecTE edges
    vecSPANDIR = fcnGLOBSTAR(repmat([0 1 0], length(vecTEDVE)), matROTANG(vecTEDVE,:)); % Spanwise direction for each HDVE (may change with rotor stuff)
end

matD = fcnDWING9(strATYPE, matEATT, matPLEX, valNELE, matELOC, matELST, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecSYM, matROTANG, vecSPANDIR, vecTEDVE);

valDLEN = length(matD);

%% Preparing to timestep
matWADJE = [];
matWELST = [];
matWDVE = [];
valWNELE = [];
matWEATT = [];
matWEIDX = [];
matWELOC = [];
matWPLEX = [];
matWDVECT = [];
matWALIGN = [];
matWVATT = [];
matWVNORM = [];
matWCENTER = [];
matWAKEGEOM = [];
valWSIZE = [];
matWCOEFF = [];
matWVLST = [];
matWROTANG = [];

% Building wing resultant
vecR = fcnRWING(valDLEN, 0, matCENTER, matDVECT, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER);

% Solving for wing coefficients
[matCOEFF] = fcnSOLVED(matD, vecR, valNELE);

for valTIMESTEP = 1:valMAXTIME
    %% Timestep to solution
    %   Move wing
    %   Generate new wake elements
    %   Create and solve WD-Matrix for new elements
    %   Solve wing D-Matrix with wake-induced velocities
    %   Solve entire WD-Matrix
    %   Relaxation procedure (Relax, create W-Matrix and W-Resultant, solve W-Matrix)
    %   Calculate surface normal forces
    %   Calculate DVE normal forces
    %   Calculate induced drag
    
    valWSIZE = length(vecTE);
    
    % Moving the wing
    [matVLST, matCENTER, matNEWWAKE] = fcnMOVEWING(matUINF, valDELTIME, matVLST, matCENTER, matELST, vecTE);
    
    % Generating new wake elements
    if any(vecTE)
        [matWAKEGEOM, matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC,...
            matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWCOEFF, matWROTANG] = fcnCREATEWAKE(valTIMESTEP, matNEWWAKE, matWAKEGEOM, matCOEFF, valWSIZE, ...
            vecTE, vecTEDVE, matPLEX, matELOC, matELST, matDVE, matWCOEFF);
        
        % Rebuild wing resultant
        vecR = fcnRWING(valDLEN, valTIMESTEP, matCENTER, matDVECT, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER);
        matCOEFF = fcnSOLVED(matD, vecR, valNELE);
        
        %         % Resolving wake D-matrix (steady)
        %         vecWLEDVE = [(valWNELE - 2*valWSIZE + 1):(valWNELE - valWSIZE)]'; % Post trailing edge row of wake HDVEs
        %         vecWLE = matWEIDX(vecWLEDVE,1);
        %         %         matNEWWAKECOEFF = fcnDWAKENEW(valWNELE, matPLEX, vecTEDVE, valWSIZE, matWPLEX, matELOC, vecTE, vecWLEDVE, matCOEFF, matWELOC, vecWLE, matDVE, matELST, matWDVE, matWELST, matWEATT, matWEIDX);
        %         %         matWCOEFF = repmat(matNEWWAKECOEFF,valTIMESTEP,1);
        
%         if flagRELAX == 1 && valTIMESTEP > 2
%             [matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, matWDVECT, matWVATT, matWVNORM, matWCENTER, matWROTANG] = ...
%                 fcnRELAX(vecUINF, valDELTIME, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, ...
%                 matROTANG, matWROTANG, matCENTER, matWCENTER);
%             
%             %                     % Resolving wake D-matrix (steady)
%             %                     vecWLEDVE = [(valWNELE - 2*valWSIZE + 1):(valWNELE - valWSIZE)]'; % Post trailing edge row of wake HDVEs
%             %                     vecWLE = matWEIDX(vecWLEDVE,1);
%             %                     matNEWWAKECOEFF = fcnDWAKENEW(valWNELE, matPLEX, vecTEDVE, valWSIZE, matWPLEX, matELOC, vecTE, vecWLEDVE, matCOEFF, matWELOC, vecWLE, matDVE, matELST, matWDVE, matWELST, matWEATT, matWEIDX);
%             %                     matWCOEFF = repmat(matNEWWAKECOEFF,valTIMESTEP,1);
%             %
%             %                     % Rebuild wing resultant
%             %                     vecR = fcnRWING(valDLEN, valTIMESTEP, matCENTER, matDVECT, matUINF, valWNELE, matWCOEFF, matWPLEX, valWSIZE, matWROTANG, matWCENTER);
%             %                     matCOEFF = fcnSOLVED(matD, vecR, valNELE);
%             
%         end
    end
    
    %% forces
%     for i = 1:valNELE
%         corners = fcnGLOBSTAR(matVLST(matDVE(i,:),:) - matCENTER(i,:), [repmat(matROTANG(i,1),3,1) repmat(matROTANG(i,2),3,1) repmat(matROTANG(i,3),3,1)]);
%         points = fcnPOLYGRID(corners(:,1), corners(:,2), 20);
%         circ = sum([0.5.*points(:,2).^2 points(:,2) 0.5.*points(:,1).^2 points(:,1) points(:,1).*points(:,2) ones(size(points(:,1)))].*matCOEFF(i,:),2);
%         
%         xi_1 = permute(matPLEX(1,1,i),[3 2 1]);
%         xi_2 = permute(matPLEX(2,1,i),[3 2 1]);
%         xi_3 = permute(matPLEX(3,1,i),[3 2 1]);
% 
%         eta_1 = permute(matPLEX(1,2,i),[3 2 1]);
%         eta_2 = permute(matPLEX(2,2,i),[3 2 1]);
%         eta_3 = permute(matPLEX(3,2,i),[3 2 1]);
%         
%         le_eta = @(x) eta_2 + (x - xi_2).*((eta_3 - eta_2)./(xi_3 - xi_2));
%         te_eta = @(x) eta_1 + (x - xi_1).*((eta_3 - eta_1)./(xi_3 - xi_1));
%         
%         f_circ = scatteredInterpolant(points(:,1), points(:,2), circ, 'linear');
% %         circ_total = integral2(@(x,y) f_circ(x,y).*-valDENSITY.*norm(matUINF(i,:)), xi_1, xi_3, te_eta, le_eta);
%         circ_total(valTIMESTEP, i) = integral2(@(x,y) f_circ(x,y).*-valDENSITY.*norm(matUINF(i,:)), xi_1, xi_3, te_eta, le_eta);
% 
%         
%         
%     end
        
end


%% Plot
[hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], matUINF, matROTANG, [3 1 4 4], 'opengl');
[hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), matUINF, matROTANG, 'r', 20);

if valTIMESTEP > 0
    [hFig1] = fcnPLOTWAKE(1, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER);
    [hFig1] = fcnPLOTCIRC(hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, matWPLEX, real(matWCOEFF), vecUINF, matWROTANG, 'b', 20);
end

view([33, 28])
toc



