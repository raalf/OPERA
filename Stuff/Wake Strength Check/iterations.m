clc
clear

hFig50 = figure(50);
clf(50);

hold on

load('QS\m5_dt1.mat');
s_t = [1:valMAXTIME].*valDELTIME.*2;
plot(s_t, vecTSITER, '-.ok');

load('QS\m5_dt0.1.mat');
s_t = [1:valMAXTIME].*valDELTIME.*2;
plot(s_t, vecTSITER, '-sk');

load('QS2\m5_dt1.mat');
s_t = [1:valMAXTIME].*valDELTIME.*2;
plot(s_t, vecTSITER, '-.^b');

load('QS2\m5_dt0.1.mat');
s_t = [1:valMAXTIME].*valDELTIME.*2;
plot(s_t, vecTSITER, '-db');

hold off

xlabel('Distanced Travelled in Semi-Chords');
ylabel('Iterations Per Timestep');

grid minor
box on
axis tight

legend('QSI (\Delta_T = 1)', 'QSI (\Delta_T = 0.1)', 'QSII (\Delta_T = 1)', 'QSII (\Delta_T = 0.1)', 'Location', 'NorthEast')

title('Iterations Per Timestep (Kussner Response, M = 5)');