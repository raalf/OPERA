clc
clear

load('chord.mat');
load('twist.mat');

hFig20 = figure(20);
clf(20);

plot(chord(:,1), chord(:,2), '-ok');
axis tight
xlim([0 1]);
ylabel('Normalized Chord Length, c/R');
xlabel('r/R');
box on
grid on

yyaxis right

plot(twist(:,1), twist(:,2), '--sb');
axis tight
ylabel('Pitch Angle, \beta (deg)');
xlim([0 1]);
legend('Chord','Pitch','Location','NorthEast');

WH = [4.5 5];
fcnFIG2LATEX(hFig20, 'tmotor_geom.pdf', WH)



