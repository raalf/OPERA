clear
clc

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
% strFILE = 'inputs/simple_wing.dat';
% strFILE = 'inputs/standard_cirrus.dat';
strFILE = 'inputs/2dve.dat';
% strFILE = 'inputs/4dve.dat';
% strFILE = 'inputs/nonplanar.dat';

[matPOINTS, strATYPE, vecSYM, flagRELAX, valMAXTIME, valDELTIME, seqALPHA, seqBETA, matTEPOINTS, matLEPOINTS] = fcnOPREAD(strFILE);

[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
    matELOC, matPLEX, matDVECT, matALIGN, matVATT, matVNORM, matCENTER, matROTANG, matVSCOMB] = fcnTRIANG(strATYPE, matPOINTS);

valTIMESTEP = 0;
flagRELAX = 0;

valMAXTIME = 0;
valALPHA = 20

[vecTE, vecLE] = fcnTELE(matTEPOINTS, matLEPOINTS, matVLST, matELST);

%% D-Matrix Creation
vecTEDVE = [];
vecSPANDIR = [];

vecLEDVE = nonzeros(sort(matEATT(vecLE,:),2,'descend'));
if ~isempty(vecTE)
    vecTEDVE = nonzeros(sort(matEATT(vecTE,:),2,'descend')); % A vector of trailing edge HDVEs, which corresponds to vecTE edges
    vecSPANDIR = fcnGLOBSTAR(repmat([0 1 0], length(vecTEDVE)), matROTANG(vecTEDVE,1), matROTANG(vecTEDVE,2), matROTANG(vecTEDVE,3)); % Spanwise direction for each HDVE (may change with rotor stuff)
end

matD = fcnDWING8(strATYPE, matEATT, matPLEX, valNELE, matELOC, matELST, matALIGN, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecLEDVE, vecSYM, matVATT, vecTEDVE, vecSPANDIR, matROTANG, matVNORM, matVSCOMB);

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

vecUINF = fcnUINFWING(valALPHA, 0);

% Building wing resultant
vecR = fcnRWING(strATYPE, valDLEN, 0, matELST, matCENTER, matDVECT, vecUINF, vecLE, vecLEDVE, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, [], matVNORM, matVLST);

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
    [matVLST, matCENTER, matNEWWAKE] = fcnMOVEWING(valALPHA, valBETA, valDELTIME, matVLST, matCENTER, matELST, vecTE);
    
    % Generating new wake elements
    if any(vecTE)
        [matWAKEGEOM, matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC,...
            matWPLEX, matWDVECT, matWALIGN, matWVATT, matWVNORM, matWCENTER, matWCOEFF, matWROTANG, matWVSCOMB] = fcnCREATEWAKE(valTIMESTEP, matNEWWAKE, matWAKEGEOM, ...
            matCOEFF, valWSIZE, matWCOEFF, vecTE, matEATT, vecTEDVE, vecSPANDIR, valWNELE, matPLEX, matELOC, matVLST, matELST, matDVE);
        
        % Rebuild wing resultant
        vecR = fcnRWING(strATYPE, valDLEN, valTIMESTEP, matEATT, matCENTER, matDVECT, vecUINF, vecTE, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, matWROTANG, matVNORM, matVLST, matWVSCOMB, matWCENTER);
        matCOEFF = fcnSOLVED(matD, vecR, valNELE);
        
        % Resolving wake D-matrix (steady)
        vecWLEDVE = [(valWNELE - 2*valWSIZE + 1):(valWNELE - valWSIZE)]'; % Post trailing edge row of wake HDVEs
        vecWLE = matWEIDX(vecWLEDVE,1);
        matNEWWAKECOEFF = fcnDWAKENEW(valWNELE, matPLEX, vecTEDVE, valWSIZE, matWPLEX, matELOC, vecTE, vecWLEDVE, vecSPANDIR, matCOEFF, matWELOC, vecWLE, matDVE, matELST, matWDVE, matWELST, matWEATT, matWCOEFF, matWALIGN, matWEIDX);
        matWCOEFF = repmat(matNEWWAKECOEFF,valTIMESTEP,1);
        
        ABC(:,:,valTIMESTEP) = matCOEFF;
    end
end

%% Plot

[hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], vecUINF, matROTANG, [3 1 4 4], 'opengl');
[hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), vecUINF, matROTANG, 'r', 10);

% granularity = .05;
% x = -0.6:granularity:1.2;
% % y = -1.2:granularity:1.2;
% y = ones(size(x)) - 0.5;
% z = -0.2:granularity:0.2;
% [X,Y,Z] = meshgrid(x,y,z);
% fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');
% 
% [s_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG, matVSCOMB, matCENTER);
% 
% % q_ind = s_ind + repmat(vecUINF, length(s_ind(:,1)),1);
% q_ind = s_ind;
% hold on
% quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3))
% hold off

%% End

if valTIMESTEP > 0
    hFig3 = figure(3);
    clf(3);
    dve = 1;
    hold on
    plot(1:valTIMESTEP, reshape(ABC(dve,1,:),1,[],1),'--xb')
    plot(1:valTIMESTEP, reshape(ABC(dve,2,:),1,[],1),'-.^b')
    
    plot(1:valTIMESTEP, reshape(ABC(dve,3,:),1,[],1),'--*r')
    plot(1:valTIMESTEP, reshape(ABC(dve,4,:),1,[],1),'-.or')
    
    plot(1:valTIMESTEP, reshape(ABC(dve,5,:),1,[],1),'->m')
    hold off
    
    box on
    grid on
    axis tight
    
    legend('A_1','A_2','B_1','B_2','C_3')
end







