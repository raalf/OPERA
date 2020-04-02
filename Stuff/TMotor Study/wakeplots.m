clc
clear

% load('Alpha 0 Results/TMotor_Relaxed_J0.3_0.0005.mat');
% 
% cd ./../../
% 
% hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, [], 'opengl');
% hFig1 = fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, valWSIZE, valPRESTEPS, matWVGRID, vecWDVESURFACE);
% view([90 0]);
% cd('G:\GIT\opera\Stuff\TMotor Study')
% 
% xlabel('X-Direction (m)');
% ylabel('Y-Direction (m)');
% axis off
% 
% WH = [4.5 3.5];
% % fcnFIG2LATEX(gcf, 'front_view_wakeplot_0', WH)


% %%
% clc
% clear
% 
% load('Alpha 0 Results/TMotor_Relaxed_J0.3_0.0005.mat');
% 
% % matVLST(:,1) = matVLST(:,1).*-1;
% % matWVLST(:,1) = matWVLST(:,1).*-1;
% % 
% % matVLST(:,2) = matVLST(:,2).*-1;
% % matWVLST(:,2) = matWVLST(:,2).*-1;
% 
% cd ./../../
% hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, [], 'opengl');
% hFig1 = fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, valWSIZE, valPRESTEPS, matWVGRID, vecWDVESURFACE);
% view([0 90])
% cd('G:\GIT\opera\Stuff\TMotor Study')
% 
% xlim([-1.9 -0.5]);
% ylim([-0.3 0.3]);
% zlim([-0.4 0.1]);
% 
% WH = [4.5*2 5];
% fcnFIG2LATEX(hFig1, 'top_down_wakeplot.pdf', WH)

%%
clc 
clear

load('Alpha 15 Results/New/TMotor_Fixed_J0.2113_0.00025_m5.mat');

cd ./../../

hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, [], 'opengl');
hFig1 = fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, valWSIZE, valPRESTEPS, matWVGRID, vecWDVESURFACE);
view([90 15]);

cd('G:\GIT\opera\Stuff\TMotor Study')

xlim([0.75 1.4]);
ylim([-0.3 0.3]);
zlim([0.15 0.35]);

xlabel('X-Direction (m)');
ylabel('Y-Direction (m)');
% axis off

WH = [4.5 3.5];
% fcnFIG2LATEX(gcf, 'front_view_wakeplot_15', WH)