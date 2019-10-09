clc
clear

fileList = dir('Alpha 90 Results/TMotor_Relaxed_*.mat');
for i = 1:size(fileList,1)
    load(['Alpha 90 Results/', fileList(i).name], 'CT', 'J', 'jj')
    CT_all(:,i) = CT;
    J_all(:,i) = J(jj);
end

J_all = [J_all(end) J_all(1:end-1)];
CT_all = [CT_all(:,end) CT_all(:,1:end-1)];

CT_all = mean(CT_all(end-10:end,:),1) + 0.0005;

hFig1 = figure(1);
clf(1);

% title('TMotor 18in, Alpha 90')

fw_fixed = [0.03183098862	0.005658531172;...
0.04774648293	0.006003732884;...
0.06366197724	0.005882723846;...
0.07957747155	0.005435878267;...
0.09549296586	0.004756671798;...
0.1114084602	0.003908846876;...
0.1273239545	0.002937200099;...
0.1432394488	0.001865059645;...
0.1591549431	0.0007197984282];

fw_relaxed = [0	0.008421754105;...
0.03501408748	0.008110173799;...
0.07002817496	0.00699369991;...
0.1050422624	0.005472895299;...
0.1400563499	0.003594867235;...
0.1750704374	0.001506334943];

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
% plot(fw_fixed(:,1), fw_fixed(:,2), '--^b')
% plot(fw_relaxed(:,1), fw_relaxed(:,2), '--sm')
plot(tunnel(:,1), tunnel(:,2), '--^r')

plot(J_all(2:end), CT_all(end,2:end), '-.sb');
grid minor
box on
axis tight
xlabel('Rotor Advance Ratio, \mu');
ylabel('Rotor Thrust Coefficient');
hold off

% hFig1 = figure(1);
% clf(1);
% plot(CT_all, '-k');
% grid minor
% box on
% axis tight
% xlabel('Timestep');
% ylabel('C_T');

CT_all = [];
J_all = [];
fileList = dir('Alpha 90 Results/TMotor_Fixed_*.mat');
for i = 1:size(fileList,1)
    load(['Alpha 90 Results/', fileList(i).name], 'CT', 'valJ')
    CT_all(:,i) = CT(end);
    J_all(:,i) = valJ;
end

hold on
plot(J_all, CT_all, '--sk');
hold off

legend('RU Test Data (3000 RPM)', 'DDE Method (Relaxed Wake)', 'Location','SouthWest','LineWidth',1)
fcnFIG2LATEX(hFig1, 'TMotor_ct.pdf', [8 5])