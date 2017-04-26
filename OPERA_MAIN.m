clear
clc

tic

% warning off

%% Header

disp('=====================================================================================');
disp('+---------------+                                      ');
disp('| RYERSON       |  OPERA V1.0  []                         []  ');
disp('| APPLIED       |              ||   ___     ___     ___   ||  ');
disp('| AERODYNAMICS  |              ||  /   \   /| |\   /   \  ||  ');
disp('| LABORATORY OF |              || |  O  |__|| ||__|  O  | ||  ');
disp('| FLIGHT        |              ||  \___/--/^^^^^\--\___/  ||  ');
disp('+---------------+              ||________|       |________||  __');
disp('        .-----------------/  \-++--------|   .   |--------++-/  \-----------------. ');
disp('       /.---------________|  |___________\__(*)__/___________|  |________---------.\');
disp('                 |    |   ''$$''   |                       |   ''$$''   |    |       ');
disp('                (o)  (o)        (o)                     (o)        (o)  (o)         ');
disp('=====================================================================================');

%% Preamble
% 
% strFILE = 'inputs/simple_wing.dat';
% strFILE = 'inputs/standard_cirrus.dat';
% strFILE = 'inputs/2dve.dat';
strFILE = 'inputs/4dve.dat';

[matPOINTS, strATYPE, vecSYM, flagRELAX, valMAXTIME, valDELTIME, seqALPHA, seqBETA, matTEPOINTS, matLEPOINTS] = fcnOPREAD(strFILE);

[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
    matELOC, matPLEX, matDVECT, matALIGN, matVATT, matVNORM, matCENTER, matROTANG] = fcnTRIANG(strATYPE, matPOINTS);

[vecTE, vecLE] = fcnTELE(matTEPOINTS, matLEPOINTS, matVLST, matELST);

% Used for STL file reading
% [TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matELOC, ...
%     matPLEX, matDVECT, matALIGN, matVATT, matVNORM, matCENTER, matROTANG] = fcnIMPORTGEOM(strFILE, strATYPE);

%% D-Matrix Creation
vecTEDVE = [];
vecSPANDIR = [];
if ~isempty(vecTE)
    vecTEDVE = nonzeros(sort(matEATT(vecTE,:),2,'descend')); % A vector of trailing edge HDVEs, which corresponds to vecTE edges
    vecSPANDIR = fcnGLOBSTAR(repmat([0 1 0], length(vecTEDVE)), matROTANG(vecTEDVE,1), matROTANG(vecTEDVE,2), matROTANG(vecTEDVE,3)); % Spanwise direction for each HDVE (may change with rotor stuff)
end

matD = fcnDWING9(strATYPE, matEATT, matPLEX, valNELE, matELOC, matELST, matALIGN, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecSYM, matVATT, vecTEDVE, vecSPANDIR, matROTANG, matVNORM);

valDLEN = length(matD);

%% Alpha Loop
for ai = 1:length(seqALPHA)
    valALPHA = deg2rad(seqALPHA(ai));
    for bi = 1:length(seqBETA)
        valBETA = deg2rad(seqBETA(bi));
        
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
        
        vecUINF = fcnUINFWING(valALPHA, valBETA);
        
        % Building wing resultant
        vecR = fcnRWING(strATYPE, valDLEN, 0, matEATT, matCENTER, matDVECT, vecUINF, vecTE, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, [], matVNORM, matVLST);
        
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
                    matWPLEX, matWDVECT, matWALIGN, matWVATT, matWVNORM, matWCENTER, matWCOEFF, matWROTANG] = fcnCREATEWAKE(valTIMESTEP, matNEWWAKE, matWAKEGEOM, ...
                    matCOEFF, valWSIZE, matWCOEFF, vecTE, matEATT, vecTEDVE, vecSPANDIR, valWNELE, matPLEX, matELOC, matVLST, matELST, matDVE);
                
                % Rebuild wing resultant
                vecR = fcnRWING(strATYPE, valDLEN, valTIMESTEP, matEATT, matCENTER, matDVECT, vecUINF, vecTE, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, matWROTANG, matVNORM, matVLST);
                matCOEFF = fcnSOLVED(matD, vecR, valNELE);
                
                % Resolving wake D-matrix (steady)
                vecWLEDVE = [(valWNELE - 2*valWSIZE + 1):(valWNELE - valWSIZE)]'; % Post trailing edge row of wake HDVEs
                vecWLE = matWEIDX(vecWLEDVE,1);
                matNEWWAKECOEFF = fcnDWAKENEW(valWNELE, matPLEX, vecTEDVE, valWSIZE, matWPLEX, matELOC, vecTE, vecWLEDVE, vecSPANDIR, matCOEFF, matWELOC, vecWLE, matDVE, matELST, matWDVE, matWELST, matWEATT, matWCOEFF, matWALIGN, matWEIDX);
                matWCOEFF = repmat(matNEWWAKECOEFF,valTIMESTEP,1);
                
                if flagRELAX == 1
                    [matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, matWDVECT, matWALIGN, matWVATT, matWVNORM, matWCENTER, matWROTANG] ...
                        = fcnRELAX(vecUINF, valDELTIME, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, matROTANG, matWROTANG);
                  
                    % Rebuild wing resultant
                    vecR = fcnRWING(strATYPE, valDLEN, valTIMESTEP, matEATT, matCENTER, matDVECT, vecUINF, vecTE, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, matWROTANG, matVNORM, matVLST);
                    matCOEFF = fcnSOLVED(matD, vecR, valNELE);                  
                    
                    % Resolving wake D-matrix (steady)
                    vecWLEDVE = [(valWNELE - 2*valWSIZE + 1):(valWNELE - valWSIZE)]'; % Post trailing edge row of wake HDVEs
                    vecWLE = matWEIDX(vecWLEDVE,1);
                    matNEWWAKECOEFF = fcnDWAKENEW(valWNELE, matPLEX, vecTEDVE, valWSIZE, matWPLEX, matELOC, vecTE, vecWLEDVE, vecSPANDIR, matCOEFF, matWELOC, vecWLE, matDVE, matELST, matWDVE, matWELST, matWEATT, matWCOEFF, matWALIGN, matWEIDX);
                    matWCOEFF = repmat(matNEWWAKECOEFF,valTIMESTEP,1);
                end
                          
                ABC(:,:,valTIMESTEP) = matCOEFF;
            end
        end
    end
end

%% Plot

[hFig1] = fcnPLOTBODY(1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, vecUINF, matROTANG, [35 1 4 4]);
[hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), vecUINF, matROTANG, 'r');
%     q_inds = fcnSDVEVEL(matCENTER, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG);
%     q_indw = fcnWINDVEL(matCENTER, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, matWROTANG);
%     q_ind = q_inds + q_indw;
%     hold on
%     quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), q_ind(:,1)+vecUINF(1), q_ind(:,2)+vecUINF(2), q_ind(:,3)+vecUINF(3), 0.25, 'g')
%     hold off
if any(vecTE) && valMAXTIME > 0
    [hFig1] = fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER);
%     [hFig1] = fcnPLOTCIRC(hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, matWPLEX, matWCOEFF, vecUINF, matWROTANG, 'b');
end
%
granularity = .1;
x = -0.6:granularity:1.2;
y = -1.2:granularity:1.2;
% y = ones(size(x)) -.5;
z = -0.2:granularity:0.2;
[X,Y,Z] = meshgrid(x,y,z);
fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');

% temp1 = [];
% temp1(:,:,1) = matVLST(matELST(:,1),:);
% temp1(:,:,2) = matVLST(matELST(:,2),:);
% edge_midpoints = mean(temp1,3);
% fpg = [matCENTER; matVLST; edge_midpoints];

% fpg = fpg + matVLST(1,:);

% fpg = [0.5 1.5 1];

% [w_ind] = fcnWINDVEL(fpg, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, matWROTANG);
[s_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG);


q_ind = s_ind + repmat(vecUINF, length(s_ind(:,1)),1);
% q_ind = s_ind;
hold on
% quiver3(fpg(:,1), fpg(:,2), fpg(:,3), w_ind(:,1), w_ind(:,2), w_ind(:,3))
quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3))
hold off
set(gcf,'Renderer','opengl');

%% End
%

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
%
% hFig4 = figure(4);
% clf(4);
% dve = 1;
% patch('Faces',[1 2 3],'Vertices',matPLEX(:,:,dve),'FaceColor','None','LineWidth',2)
%
% axis equal
% grid on
% box on
% axis tight






