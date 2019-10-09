clc
clear

%%
hFig22 = figure(22);
clf(22);
load('Stuff/Gust Response/Kussner.mat')

subplot(2,2,1)

plot(Kussner(:,1), smooth(Kussner(:,2)), '-k', 'LineWidth', 1.5);
box on
axis tight
grid minor

hold on
load('Kussner_m1_0.5_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME');
plot(valDELTIME.*[0:valMAXTIME-1]./0.5, CL2D, '-b', 'LineWidth', 1.5);

load('Kussner_m3_0.5_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME');
plot(valDELTIME.*[0:valMAXTIME-1]./0.5, CL2D, '-.r', 'LineWidth', 1.5);

load('Kussner_m5_0.5_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME');
plot(valDELTIME.*[0:valMAXTIME-1]./0.5, CL2D, '--m', 'LineWidth', 1.5);

hold off

legend('Kussner Function', 'M = 1, $\Delta x_w/\Delta x_c = 0.25$', 'M = 3, $\Delta x_w/\Delta x_c = 0.75$', 'M = 5, $\Delta x_w/\Delta x_c = 1.25$', 'Location', 'SouthEast','Interpreter','Latex','FontSize',14)
xlabel('Distance Traversed by Gust in Semi-Chords')
ylabel('c_l')
title('\Delta_T = 0.25s, \Deltax_w/c = 0.25')

%%
% hFig23 = figure(23);
% clf(23);
% load('Stuff/Gust Response/Kussner.mat')

subplot(2,2,2)

plot(Kussner(:,1), smooth(Kussner(:,2)), '-k', 'LineWidth', 1.5);
box on
axis tight
grid minor

hold on
load('Kussner_m1_0.5_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME');
plot(valDELTIME.*[0:valMAXTIME-1]./0.5, CL2D, '-b', 'LineWidth', 1.5);

load('Kussner_m3_0.5_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME');
plot(valDELTIME.*[0:valMAXTIME-1]./0.5, CL2D, '-.r', 'LineWidth', 1.5);

load('Kussner_m5_0.5_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME');
plot(valDELTIME.*[0:valMAXTIME-1]./0.5, CL2D, '--m', 'LineWidth', 1.5);

hold off

legend('Kussner Function', 'M = 1, $\Delta x_w/\Delta x_c = 0.5$', 'M = 3, $\Delta x_w/\Delta x_c = 1.5$', 'M = 5, $\Delta x_w/\Delta x_c = 2.5$', 'Location', 'SouthEast','Interpreter','Latex','FontSize',14)
xlabel('Distance Traversed by Gust in Semi-Chords')
ylabel('c_l')
title('\Delta_T = 0.5s, \Deltax_w/c = 0.5')

%%
% hFig24 = figure(24);
% clf(24);
% load('Stuff/Gust Response/Kussner.mat')

subplot(2,2,3)

plot(Kussner(:,1), smooth(Kussner(:,2)), '-k', 'LineWidth', 1.5);
box on
axis tight
grid minor

hold on
load('Kussner_m1_1_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME');
plot(valDELTIME.*[0:valMAXTIME-1]./0.5, CL2D, '-b', 'LineWidth', 1.5);

load('Kussner_m3_1_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME');
plot(valDELTIME.*[0:valMAXTIME-1]./0.5, CL2D, '-.r', 'LineWidth', 1.5);

load('Kussner_m5_1_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME');
plot(valDELTIME.*[0:valMAXTIME-1]./0.5, CL2D, '--m', 'LineWidth', 1.5);

hold off

legend('Kussner Function', 'M = 1, $\Delta x_w/\Delta x_c = 1$', 'M = 3, $\Delta x_w/\Delta x_c = 3$', 'M = 5, $\Delta x_w/\Delta x_c = 5$', 'Location', 'SouthEast','Interpreter','Latex','FontSize',14)
xlabel('Distance Traversed by Gust in Semi-Chords')
ylabel('c_l')
title('\Delta_T = 1s, \Deltax_w/c = 1')

%%
% hFig25 = figure(25);
% clf(25);
% load('Stuff/Gust Response/Kussner.mat')

subplot(2,2,4)

plot(Kussner(:,1), smooth(Kussner(:,2)), '-k', 'LineWidth', 1.5);
box on
axis tight
grid minor

hold on
load('Kussner_m1_1.5_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME');
plot(valDELTIME.*[0:valMAXTIME-1]./0.5, CL2D, '-b', 'LineWidth', 1.5);

load('Kussner_m3_1.5_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME');
plot(valDELTIME.*[0:valMAXTIME-1]./0.5, CL2D, '-.r', 'LineWidth', 1.5);

load('Kussner_m5_1.5_fixed.mat', 'CL2D', 'valDELTIME', 'valMAXTIME');
plot(valDELTIME.*[0:valMAXTIME-1]./0.5, CL2D, '--m', 'LineWidth', 1.5);

hold off

legend('Kussner Function', 'M = 1, $\Delta x_w/\Delta x_c = 1.5$', 'M = 3, $\Delta x_w/\Delta x_c = 4.5$', 'M = 5, $\Delta x_w/\Delta x_c = 7.5$', 'Location', 'SouthEast','Interpreter','Latex','FontSize',14)
xlabel('Distance Traversed by Gust in Semi-Chords')
ylabel('c_l')
title('\Delta_T = 1.5s, \Deltax_w/c = 1.5')