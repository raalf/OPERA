clc
clear

load('Alpha 0 Results/TMotor_Relaxed_J0.3_0.0005.mat');

cd ./../../

hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, [], 'opengl');
hFig1 = fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, valWSIZE, valPRESTEPS, matWVGRID, vecWDVESURFACE);
view([90 0]);

cd('G:\GIT\opera\Stuff\TMotor Study')

xlim([0.5 1.9]);
ylim([-0.3 0.3]);
zlim([-0.4 0.1]);

xlabel('X-Direction (m)');
ylabel('Y-Direction (m)');
axis off

WH = [4.5 3.5];
% fcnFIG2LATEX(gcf, 'front_view_wakeplot_0', WH)

%%
clc 
clear

load('Alpha 15 Results/TMotor_Relaxed_J0.1115_0.0005_2.mat');

cd ./../../

hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, [], 'opengl');
hFig1 = fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, valWSIZE, valPRESTEPS, matWVGRID, vecWDVESURFACE);
view([90 15]);

cd('G:\GIT\opera\Stuff\TMotor Study')

xlim([0.35 0.8]);
ylim([-0.3 0.3]);
zlim([-0.2 0.2]);

xlabel('X-Direction (m)');
ylabel('Y-Direction (m)');
axis off

WH = [4.5 3.5];
% fcnFIG2LATEX(gcf, 'front_view_wakeplot_15', WH)