clc
clear


% hFig24 = figure(24);
% clf(24);
% load('Stuff/Gust Response/ZAERO.mat')
% plot(smooth(ZAERO(:,1)), smooth(ZAERO(:,2)), '-k');
% hold on
% 
% load('Stuff/Gust Response/UVLM.mat')
% plot(smooth(UVLM(:,1)), smooth(UVLM(:,2)), '-.k');
% 
% plot(0.135+((valDELTIME.*([0:valMAXTIME-1]))), CL2, '--b');
% box on
% axis tight
% grid minor
% %profile viewer
% 
% GOLAND_X = 0.135+((valDELTIME.*([0:valMAXTIME-1])));


hFig2 = figure(2);
clf(2);

load('ZAERO.mat')
plot(smooth(ZAERO(:,1)), smooth(ZAERO(:,2)), '-ok');
hold on
% load('UVLM.mat')
% plot(smooth(UVLM(:,1)), smooth(UVLM(:,2)), '-.k');

load('Goland_m3_fixed.mat','CL2','GOLAND_X');
plot(GOLAND_X+0.01, CL2, '--db');

% load('Goland_m4_fixed.mat','CL2','GOLAND_X');
% plot(GOLAND_X, CL2, '-sr');

load('Goland_m5_fixed.mat','CL2','GOLAND_X');
plot(GOLAND_X, CL2, '--^m');

% legend('ZAERO','M = 3', 'M = 5', 'Location','NorthEast')
legend('ZAERO', 'M = 3, $\Delta x_w/\Delta x_c = 0.8202$', 'M = 5, $\Delta x_w/\Delta x_c = 1.3670$', 'Location', 'NorthEast','Interpreter','Latex','FontSize',14)

grid minor
box on
axis tight
xlabel('Time (s)');
ylabel('C_L');
title('Goland Wing Comparison (\Deltax_w/c = 0.2734)');