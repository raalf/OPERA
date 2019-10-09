clc
clear

%%
hFig22 = figure(22);
clf(22);
load('Stuff/Gust Response/Kussner.mat')

load('Kussner_m3_0.25_fixed.mat', 'valDELTIME');

plot(Kussner(:,1).*2, smooth(Kussner(:,2)), '-k');
box on
axis tight
grid minor

hold on
load('Kussner_m2_0.25_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME', 'valGUSTSTART');
plot([1:valMAXTIME-1] - valGUSTSTART, CL2D(2:end), '-sb');

load('Kussner_m4_0.25_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME', 'valGUSTSTART');
plot([1:valMAXTIME-1] - valGUSTSTART, CL2D(2:end), '-.^r');

load('Kussner_m7_0.25_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME');
plot([1:valMAXTIME-1] - valGUSTSTART, CL2D(2:end), '--om');

hold off

legend('Kussner Function', 'M = 2', 'M = 4', 'M = 7', 'Location', 'SouthEast')
xlabel('Timestep')
ylabel('Two-Dimensional Lift Coefficient')

fcnFIG2LATEX(hFig22, 'Kussner.pdf', [8 5])