clear
clc

% %% Preamble
% %
% strFILE = 'inputs/Stock_Test1.dat';
% 
% [matPOINTS, strATYPE, vecSYM, flagRELAX, valMAXTIME, valDELTIME, seqALPHA, seqBETA, matTEPOINTS, matLEPOINTS] = fcnOPREAD(strFILE);
% 
% [TR, matADJE, matELST, matVLST, matDVE, valNELE, matEATT, matEIDX, ...
%     matELOC, matPLEX, matDVECT, matALIGN, matVATT, matVNORM, matCENTER, matROTANG, matVSCOMB] = fcnTRIANG(matPOINTS);
% 
% flagRELAX = 0;
% 
% valMAXTIME = 5;
% valDELTIME = 0.05;
% valALPHA = 0;
% 
% %% D-Matrix Creation
% vecTEDVE = [];
% vecLEDVE = [];
% vecSPANDIR = [];
% vecTE = [];
% vecLE = [];
% 
% matD = fcnDWING9(strATYPE, matEATT, matPLEX, valNELE, matELOC, matELST, matALIGN, matVLST, matCENTER, matDVE, matDVECT, vecTE, vecLE, vecLEDVE, vecSYM, matVATT, vecTEDVE, vecSPANDIR, matROTANG, matVNORM, matVSCOMB);
% 
% valDLEN = length(matD);
% 
% %% Preparing to timestep
% valTIMESTEP = 0;
% matWADJE = [];
% matWELST = [];
% matWDVE = [];
% valWNELE = [];
% matWEATT = [];
% matWEIDX = [];
% matWELOC = [];
% matWPLEX = [];
% matWDVECT = [];
% matWALIGN = [];
% matWVATT = [];
% matWVNORM = [];
% matWCENTER = [];
% matWAKEGEOM = [];
% valWSIZE = [];
% matWCOEFF = [];
% matWVLST = [];
% 
% % Stock's velocity field
% matUINF_VLST = [matVLST(:,1).*0, cos(2.*pi.*matVLST(:,2)).*cos(2.*pi.*matVLST(:,3)), sin(2.*pi.*matVLST(:,2)).*sin(2.*pi.*matVLST(:,3))];
% matUINF = [matCENTER(:,1).*0, cos(2.*pi.*matCENTER(:,2)).*cos(2.*pi.*matCENTER(:,3)), sin(2.*pi.*matCENTER(:,2)).*sin(2.*pi.*matCENTER(:,3))];
% 
% % Building wing resultant
% vecR = fcnRWING(strATYPE, valDLEN, 0, matELST, matCENTER, matDVECT, matUINF, vecLE, vecLEDVE, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, [], matVNORM, matVLST);
% 
% % Solving for wing coefficients
% [matCOEFF] = fcnSOLVED(matD, vecR, valNELE);
% 
% % matCOEFF = repmat([0 1 0 0 0 0], valNELE, 1);
% 
% % vidx = matVLST(:,2) > 0 & matVLST(:,2) < 1;
% vidx = true(size(matVLST,1),1);
% q_ind_VLST = matUINF_VLST.*0;
% q_ind = matCENTER.*0;
% for valTIMESTEP = 1:valMAXTIME
% 
% matVLST(vidx,:) = matVLST(vidx,:) + (matUINF_VLST(vidx,:) + q_ind_VLST(vidx,:)).*valDELTIME;
% 
% matUINF_VLST = [matVLST(:,1).*0, cos(2.*pi.*matVLST(:,2)).*cos(2.*pi.*matVLST(:,3)), sin(2.*pi.*matVLST(:,2)).*sin(2.*pi.*matVLST(:,3))];
% matUINF = [matCENTER(:,1).*0, cos(2.*pi.*matCENTER(:,2)).*cos(2.*pi.*matCENTER(:,3)), sin(2.*pi.*matCENTER(:,2)).*sin(2.*pi.*matCENTER(:,3))];
% 
% matGEOM = [];
% matGEOM(:,:,1) = matVLST(matDVE(:,1,1),:);
% matGEOM(:,:,2) = matVLST(matDVE(:,2,1),:);
% matGEOM(:,:,3) = matVLST(matDVE(:,3,1),:);
% 
% [~, ~, matELST, ~, ~, ~, ~, matEIDX, matELOC, matPLEX, ...
%     matDVECT, matALIGN, matVATT, matVNORM, matCENTER, matROTANG] = fcnTRIANG(matGEOM);
% 
% vecR = fcnRWING(strATYPE, valDLEN, 0, matELST, matCENTER, matDVECT, matUINF, vecLE, vecLEDVE, valWNELE, matWCOEFF, matWDVE, matWDVECT, matWVLST, matWPLEX, valWSIZE, [], matVNORM, matVLST);
% [matCOEFF] = fcnSOLVED(matD, vecR, valNELE);
% 
% q_ind_VLST = fcnSDVEVEL(matVLST, valNELE, matCOEFF, matDVE, matDVECT, matVLST, matPLEX, matROTANG, matVSCOMB, matCENTER);
% 
% end
% 
% %% Plot
% 
% [hFig1] = fcnPLOTBODY(1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, [], matUINF, matROTANG, [3 1 4 4], 'opengl');
% % [hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, real(matCOEFF), matUINF, matROTANG, 'r', 5);
% 
% hold on
% % quiver3(matVLST(:,1), matVLST(:,2), matVLST(:,3), matUINF_VLST(:,1), matUINF_VLST(:,2), matUINF_VLST(:,3))
% 
% h1 = quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), matDVECT(:,1,1), matDVECT(:,2,1), matDVECT(:,3,1), 0.25, 'k'); % xsi
% h2 = quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), matDVECT(:,1,2), matDVECT(:,2,2), matDVECT(:,3,2), 0.25, 'b'); % eta
% h3 = quiver3(matCENTER(:,1),matCENTER(:,2),matCENTER(:,3), matDVECT(:,1,3), matDVECT(:,2,3), matDVECT(:,3,3), 0.25, 'm'); % zeta (normal)
% 
% hold off
% 

%%
clc
clear

load('temp1.mat')


v_row = 181:200;

ii = 1;
iii = 1;
for jj = 1:length(v_row)
    
    [~, b] = ismembertol(v_row(jj), matDVE,'OutputAllIndices',true);

    [elem(:,1) elem(:,2)] = ind2sub(size(matDVE), cell2mat(b));

    for i = 1:size(elem,1)
        point = matPLEX(elem(i,2),:,elem(i,1));
        vort = [2.*matCOEFF(elem(i,1),3).*point(:,1) + matCOEFF(elem(i,1),4) + matCOEFF(elem(i,1),5).*point(:,2), 2.*matCOEFF(elem(i,1),1).*point(:,2) + matCOEFF(elem(i,1),2) + matCOEFF(elem(i,1),5).*point(:,1), point(:,2).*0];
        vort_glob(ii,:) = fcnSTARGLOB(vort, matROTANG(elem(i,1),1), matROTANG(elem(i,1),2), matROTANG(elem(i,1),3));
        point_glob(ii,:) = matVLST(v_row(jj),:);
        ii = ii + 1;
    end
    
    vort_avg(iii,:) = mean(vort_glob(ii-size(elem,1):ii-1,:),1);
    point_avg(iii,:) = point_glob(ii-1,:);
    
    iii = iii + 1;
    elem = [];
    point = [];
    vort = [];
end

hFig28 = figure(28);
clf(28);

scatter(point_glob(:,2), vort_glob(:,1));
% scatter(point_avg(:,2), vort_avg(:,1))
grid minor
box on
axis tight

ylabel('Gamma_x','FontSize',15);
xlabel('Y-Location (m)','FontSize',15);









