clc
clear

hFig23 = figure(23);
clf(23);
load('Stuff/Gust Response/Kussner.mat')
plot(Kussner(:,1), smooth(Kussner(:,2)), '-k', 'LineWidth', 1.5);
box on
axis tight
grid minor

hold on
load('Kussner_m3_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME');
plot(valDELTIME.*[0:valMAXTIME-1]./0.5, CL2D, '-b', 'LineWidth', 1.5);

load('Kussner_m4_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME');
plot(valDELTIME.*[0:valMAXTIME-1]./0.5, CL2D, '-.r', 'LineWidth', 1.5);

load('Kussner_m6_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME');
plot(valDELTIME.*[0:valMAXTIME-1]./0.5, CL2D, '--m', 'LineWidth', 1.5);

hold off

legend('Kussner Function', 'M = 3', 'M = 4', 'M = 6', 'Location', 'SouthEast')
xlabel('Distance Traversed by Gust in Semi-Chords')
ylabel('C_l')
title('Unsteady DDE Response to Sharp-Edged Gust')