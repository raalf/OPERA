clc
clear

valRPM = 3000;

offset = 0;
CT_all = [];
J_all = [];
fileList = dir('Alpha 0 Results/TMotor_Relaxed_*.mat');
for i = 1:size(fileList,1)
    load(['Alpha 0 Results/', fileList(i).name], 'CT', 'valJ', 'valDELTIME')
    num_azm = 1/((valRPM/60).*valDELTIME);
    CT = CT(~isnan(CT));
    CT_all(:,i) = mean(CT(end-(num_azm*2):end));
    J_all(:,i) = valJ;
end

hFig2 = figure(2);
clf(2);

tunnel = [0.0416	0.0098;...
0.0593	0.0101;...
0.0766	0.0106;...
0.0935	0.0111;...
0.1113	0.0117;...
0.1298	0.0124;...
0.1496	0.013;...
0.1695	0.0135;...
0.1873	0.014;...
0.2055	0.0144;...
0.2231	0.0147;...
0.2405	0.0151;...
0.2579	0.0155;...
0.2753	0.0159;...
0.2915	0.0163;...
0.3087	0.0166;...
0.3249	0.0172];

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
fileList = dir('Alpha 0 Results/TMotor_Fixed_*.mat');
for i = 1:size(fileList,1)
    load(['Alpha 0 Results/', fileList(i).name], 'CT', 'valJ', valDELTIME)
    num_azm = 1/((valRPM/60).*valDELTIME);
    CT = CT(~isnan(CT));
    CT_all(:,i) = mean(CT(end-(num_azm*2):end));
    J_all(:,i) = valJ;
end

hold on
plot(J_all, CT_all + offset, '--sk');
hold off

legend('RU Test Data (3000 RPM)', 'DDE Method (Relaxed Wake)', 'Location','NorthWest','LineWidth',1)
fcnFIG2LATEX(hFig2, 'TMotor_ct_0.pdf', [8 5])

