clc
clear

%%
hFig22 = figure(22);
clf(22);
load('Stuff/Gust Response/Kussner.mat')

load('Kussner_m2_0.5_fixed.mat', 'valDELTIME');

plot(Kussner(:,1), Kussner(:,2), '-k');
box on
axis tight
grid minor

hold on
load('Kussner_m2_0.5_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME', 'valGUSTSTART');
plot(([1:valMAXTIME-1] - valGUSTSTART).*valDELTIME.*2, CL2D(2:end), '-sb');

load('Kussner_m4_0.25_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME', 'valGUSTSTART');
plot(([1:valMAXTIME-1] - valGUSTSTART).*valDELTIME.*2, CL2D(2:end), '-.^r');

load('Kussner_m8_0.125_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME');
plot(([1:valMAXTIME-1] - valGUSTSTART).*valDELTIME.*2, CL2D(2:end), '--m');

hold off

legend('Kussner Function', 'M = 2', 'M = 4', 'M = 8', 'Location', 'SouthEast')
xlabel('Distance Travelled by Gust in Semi-Chords')
ylabel('Two-Dimensional Lift Coefficient')

fcnFIG2LATEX(hFig22, 'Kussner.pdf', [8 5])