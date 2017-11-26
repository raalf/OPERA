clear
clc

tic

% warning off

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

% matPOINTS(:,:,1) = [1 1 0; 1 1 0];
% matPOINTS(:,:,2) = [0 0 0; 1 0 0];
% matPOINTS(:,:,3) = [0 1 0; 0 0 0];

[TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
    matELOC, matPLEX, matDVECT, matALIGN, matVATT, matVNORM, matCENTER, matROTANG, matVSCOMB] = fcnTRIANG(strATYPE, matPOINTS);

valMAXTIME = 0
valTIMESTEP = 0
seqALPHA = 10
flagRELAX = 0

vecM = 1
vecN = 10

[vecTE, vecLE] = fcnTELE(matTEPOINTS, matLEPOINTS, matVLST, matELST);

% Used for STL file reading
% strFILE = 'CAD Geom/cube.stl'
% strATYPE = 'PC'
% [TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, matELOC, ...
%     matPLEX, matDVECT, matALIGN, matVATT, matVNORM, matCENTER, matROTANG] = fcnIMPORTGEOM(strFILE, strATYPE);
% vecTE = []

%% D-Matrix Creation
vecTEDVE = [];
vecSPANDIR = [];
if ~isempty(vecTE)
    vecTEDVE = nonzeros(sort(matEATT(vecTE,:),2,'descend')); % A vector of trailing edge HDVEs, which corresponds to vecTE edges
    vecSPANDIR = fcnGLOBSTAR(repmat([0 1 0], length(vecTEDVE)), matROTANG(vecTEDVE,1), matROTANG(vecTEDVE,2), matROTANG(vecTEDVE,3)); % Spanwise direction for each HDVE (may change with rotor stuff)
end

matD = fcnDWING8(strATYPE, matEATT, matPLEX, valNELE, matELOC, matELST, matALIGN, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecSYM, matVATT, vecTEDVE, vecSPANDIR, matROTANG, matVNORM, matVSCOMB);

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
                
                if flagRELAX == 1
                    [matWADJE, matWELST, matWVLST, matWDVE, valWNELE, matWEATT, matWEIDX, matWELOC, matWPLEX, matWDVECT, matWALIGN, matWVATT, matWVNORM, matWCENTER, matWROTANG] ...
                        = fcnRELAX(vecUINF, valDELTIME, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, matROTANG, matWROTANG, matVSCOMB, matCENTER, matWVSCOMB, matWCENTER);
                  
                    % Rebuild wing resultant
                    vecR = fcnRWING(strATYPE, valDLEN, valTIMESTEP, matEATT, matCENTER, matDVECT, vecUINF, vecTE, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, matWROTANG, matVNORM, matVLST, matWVSCOMB, matWCENTER);
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

[hFig1] = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], vecUINF, matROTANG, [3 1 4 4], 'opengl');

% dve_num = 1;
% dir = 1;
% for i = 1:3  
%     edge_points = matVLST(matELST(matEIDX(dve_num,i),:),:);
%     hold on
%     if matVSCOMB(dve_num,i,dir) == 1
%         plot3(edge_points(:,1), edge_points(:,2), edge_points(:,3),'--r','LineWidth',2)
%     elseif matVSCOMB(dve_num,i,dir) == -1
%         plot3(edge_points(:,1), edge_points(:,2), edge_points(:,3),'--b','LineWidth',2)
%     end
%     hold off
% end

% % Plotting global and local to visualize
% test_num = 2;
% hFig2 = figure(2);
% clf(2);
% subplot(2,1,1)
% patch(matPLEX(:,1,test_num),matPLEX(:,2,test_num),'b')
% alpha(0.5);
% xlabel('u-direction','FontSize',15);
% ylabel('v-direction','FontSize',15);
% axis equal
% grid minor
% box on
% subplot(2,1,2)
% P = [matVLST(matDVE(test_num,1),:); matVLST(matDVE(test_num,2),:); matVLST(matDVE(test_num,3),:)];
% patch(P(:,1),P(:,2),P(:,3),'r')
% alpha(0.5);
% xlabel('X-direction','FontSize',15);
% ylabel('Y-direction','FontSize',15);
% zlabel('Z-direction','FontSize',15);
% axis equal
% grid minor
% box on

[hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), vecUINF, matROTANG, 'r');
% %     q_inds = fcnSDVEVEL(matCENTER, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG);
% %     q_indw = fcnWINDVEL(matCENTER, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, matWROTANG);
% %     q_ind = q_inds + q_indw;
% %     hold on
% %     quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), q_ind(:,1)+vecUINF(1), q_ind(:,2)+vecUINF(2), q_ind(:,3)+vecUINF(3), 0.25, 'g')
% %     hold off
if any(vecTE) && valMAXTIME > 0
    [hFig1] = fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER);
%     [hFig1] = fcnPLOTCIRC(hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, matWPLEX, matWCOEFF, vecUINF, matWROTANG, 'b');
end

% % %
% granularity = .05;
% x = -0.6:granularity:1.2;
% % y = -1.2:granularity:1.2;
% y = ones(size(x)) - 0.5;
% z = -0.2:granularity:0.2;
% [X,Y,Z] = meshgrid(x,y,z);
% fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');
% 
% % temp1 = [];
% % temp1(:,:,1) = matVLST(matELST(:,1),:);
% % temp1(:,:,2) = matVLST(matELST(:,2),:);
% % edge_midpoints = mean(temp1,3);
% % fpg = [matCENTER; matVLST; edge_midpoints];
% % 
% % fpg = fpg + matCENTER;
% 
% % fpg = [0.5 0 1];
% 
% % [w_ind] = fcnWINDVEL(fpg, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, matWROTANG, matVSCOMB);
% [s_ind] = fcnSDVEVEL(fpg, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG, matVSCOMB, matCENTER);
% 
% 
% q_ind = s_ind + repmat(vecUINF, length(s_ind(:,1)),1);
% % q_ind = s_ind;
% hold on
% % quiver3(fpg(:,1), fpg(:,2), fpg(:,3), w_ind(:,1), w_ind(:,2), w_ind(:,3))
% quiver3(fpg(:,1), fpg(:,2), fpg(:,3), q_ind(:,1), q_ind(:,2), q_ind(:,3))
% hold off
% set(gcf,'Renderer','opengl');

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






