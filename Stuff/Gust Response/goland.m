clc
clear

hFig2 = figure(2);
clf(2);

load('ZAERO.mat')
plot(smooth(ZAERO(:,1)), smooth(ZAERO(:,2)), '-k');
hold on
% load('UVLM.mat')
% plot(smooth(UVLM(:,1)), smooth(UVLM(:,2)), '-.k');
offset = 0.05
load('Goland_m2_0.025.mat','CL2','GOLAND_X');
plot(GOLAND_X+offset, CL2, '-sb');

% load('Goland_m4_fixed.mat','CL2','GOLAND_X');
% plot(GOLAND_X, CL2, '-sr');

load('Goland_m4_0.025.mat','CL2','GOLAND_X');
plot(GOLAND_X+offset, CL2, '-.^r');

load('Goland_m7_0.025.mat','CL2','GOLAND_X');
plot(GOLAND_X+offset, CL2, '--om');

legend('ZAERO', 'M = 2', 'M = 4', 'M = 7', 'Location', 'NorthEast')

grid minor
box on
axis tight
xlabel('Time (s)');
ylabel('Lift Coefficient');

fcnFIG2LATEX(hFig2, 'Goland.pdf', [8 5])


