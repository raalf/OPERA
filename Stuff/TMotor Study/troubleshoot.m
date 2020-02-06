clc
clear


hFig1 = figure(1);
clf(1);

load('Alpha 0 Results/TMotor_Fixed_J0.3_0.0005_2.mat')

plot(CT_U, '-k');

load('Alpha 15 Results/TMotor_Fixed_J0.1109_0.0005.mat')
yyaxis right
plot(CT_U, '--r');

box on
grid minor
axis tight