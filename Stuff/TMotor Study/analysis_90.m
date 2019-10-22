clc
clear

offset = 0;
CT_all = [];
J_all = [];
fileList = dir('Alpha 90 Results/TMotor_Relaxed_*.mat');
for i = 1:size(fileList,1)
    load(['Alpha 90 Results/', fileList(i).name], 'CT', 'valJ')
    CT = CT(~isnan(CT));
    CT_all(:,i) = CT(end);
    J_all(:,i) = valJ;
end

hFig1 = figure(1);
clf(1);

tunnel = [0.0519	0.0078;...
0.0582	0.0072;...
0.0639	0.0068;...
0.0704	0.0064;...
0.0766	0.0059;...
0.0827	0.0055;...
0.0892	0.0049;...
0.0957	0.0045;...
0.1026	0.0038;...
0.1093	0.0033;...
0.116	0.0027;...
0.1228	0.0021;...
0.1303	0.0014;...
0.1368	0.0007];

hold on
plot(tunnel(:,1), tunnel(:,2), '--^r')

plot(J_all, CT_all + offset, '-.sb');
grid minor
box on
axis tight
xlabel('Rotor Advance Ratio, \mu');
ylabel('Rotor Thrust Coefficient');
hold off

CT_all = [];
J_all = [];
fileList = dir('Alpha 90 Results/TMotor_Fixed_*.mat');
for i = 1:size(fileList,1)
    load(['Alpha 90 Results/', fileList(i).name], 'CT', 'valJ')
    CT_all(:,i) = CT(end);
    J_all(:,i) = valJ;
end

hold on
plot(J_all, CT_all + offset, '--sk');
hold off

legend('RU Test Data (3000 RPM)', 'DDE Method (Relaxed Wake)', 'Location','SouthWest','LineWidth',1)
fcnFIG2LATEX(hFig1, 'TMotor_ct_90.pdf', [8 5])

