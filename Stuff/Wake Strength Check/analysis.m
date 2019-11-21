clc
clear

hFig(34) = figure(34);
clf(34)
load('m5_dt1_wfix2.mat', 'matWVLST', 'vecWVMU', 'matWDVE', 'valWNELE', 'matWVLST', 'matWELST', 'matWDVECT', 'matWCENTER', 'matWPLEX', 'matWCOEFF', 'matWROTANG', 'tmp_wind', 'valDELTIME', 'valMAXTIME');
hold on
tmp_wind1 = [valDELTIME.*[1:valMAXTIME]' tmp_wind(:,3)];
scatter3(matWVLST(:,1), matWVLST(:,2), vecWVMU, 'ok')
scatter3(matWCENTER(:,1), matWCENTER(:,2), matWCOEFF(:,6), 'ok')
% fcnPLOTCIRC(hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, matWPLEX, matWCOEFF, [], matWROTANG, [], 30)
hold off

load('m5_dt0.5_wfix2.mat', 'matWVLST', 'vecWVMU', 'matWDVE', 'valWNELE', 'matWVLST', 'matWELST', 'matWDVECT', 'matWCENTER', 'matWPLEX', 'matWCOEFF', 'matWROTANG', 'tmp_wind', 'valDELTIME', 'valMAXTIME');
hold on
tmp_wind2 = [valDELTIME.*[1:valMAXTIME]' tmp_wind(:,3)];
scatter3(matWVLST(:,1), matWVLST(:,2), vecWVMU, 'sb')
scatter3(matWCENTER(:,1), matWCENTER(:,2), matWCOEFF(:,6), 'sb')
% fcnPLOTCIRC(hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, matWPLEX, matWCOEFF, [], matWROTANG, [], 30)
hold off

load('m5_dt0.25_wfix2.mat', 'matWVLST', 'vecWVMU', 'matWDVE', 'valWNELE', 'matWVLST', 'matWELST', 'matWDVECT', 'matWCENTER', 'matWPLEX', 'matWCOEFF', 'matWROTANG', 'tmp_wind', 'valDELTIME', 'valMAXTIME');
hold on
tmp_wind3 = [valDELTIME.*[1:valMAXTIME]' tmp_wind(:,3)];
scatter3(matWVLST(:,1), matWVLST(:,2), vecWVMU, 'rd')
scatter3(matWCENTER(:,1), matWCENTER(:,2), matWCOEFF(:,6), 'rd')
% fcnPLOTCIRC(hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, matWPLEX, matWCOEFF, [], matWROTANG, [], 30)
hold off

load('m5_dt0.05_wfix2.mat', 'matWVLST', 'vecWVMU', 'matWDVE', 'valWNELE', 'matWVLST', 'matWELST', 'matWDVECT', 'matWCENTER', 'matWPLEX', 'matWCOEFF', 'matWROTANG', 'tmp_wind', 'valDELTIME', 'valMAXTIME');
hold on
tmp_wind4 = [valDELTIME.*[1:valMAXTIME]' tmp_wind(:,3)];
scatter3(matWVLST(:,1), matWVLST(:,2), vecWVMU, '*m')
scatter3(matWCENTER(:,1), matWCENTER(:,2), matWCOEFF(:,6), '*m')
% fcnPLOTCIRC(hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, matWPLEX, matWCOEFF, [], matWROTANG, [], 30)
hold off

grid minor
box on
axis tight

xlabel('X-Dir');
ylabel('Y-Dir');
zlabel('Z-Dir');

legend('\Delta_T = 1s', '\Delta_T = 0.5s', '\Delta_T = 0.25s', '\Delta_T = 0.05s')


hFig5 = figure(5);
plot(tmp_wind1(:,1), tmp_wind1(:,2), '-ok');
hold on
plot(tmp_wind2(:,1), tmp_wind2(:,2), '--sb');
plot(tmp_wind3(:,1), tmp_wind3(:,2), '-.rd');
plot(tmp_wind4(:,1), tmp_wind4(:,2), '--*m');
hold off
grid minor
box on
axis tight
xlabel('Time (s)');
ylabel('Mean wake-induced z-velocity')
legend(['\Delta_T = 1s, \Deltax_w/\Deltax_c = ', num2str(1/(1/5))], ['\Delta_T = 0.5s, \Deltax_w/\Deltax_c = ', num2str(0.5/(1/5))], ['\Delta_T = 0.25s, \Deltax_w/\Deltax_c = ', num2str(0.25/(1/5))], ['\Delta_T = 0.05s, \Deltax_w/\Deltax_c = ', num2str(0.05/(1/5))],'Location','SouthEast')



