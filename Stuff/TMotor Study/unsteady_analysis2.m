clc
clear

% load('G:\GIT\opera\Stuff\TMotor Study\Tunnel Testing\2019-10-31\31-Oct-2019 12.41.35_Scorpion_ASI_T-Motor 18in_RPM3000_Alpha30_21.5736.mat') % 0.3004
load('G:\GIT\opera\Stuff\TMotor Study\Tunnel Testing\2019-11-06\06-Nov-2019 12.34.15_Scorpion_ASI_T-Motor 18in_RPM3000_Alpha0_9.667.mat') % 0.1346

close all

%% Getting OPERA data
% load('Alpha 0 Results/TMotor_Relaxed_J0.3.mat', 'CT_U', 'CT', 'valDELTIME', 'matDVE', 'matVLST', 'vecHUB', 'vecDVESURFACE', 'matDGAMMADT', 'matINTCIRC', 'vecTEDVE', 'tmpDVETHRUST', 'matSPANDIR', 'valRPM')
load('Alpha 0 Results/TMotor_Relaxed_J0.1346.mat', 'CT_U', 'CT', 'valDELTIME', 'matDVE', 'matVLST', 'vecHUB', 'vecDVESURFACE', 'matDGAMMADT', 'matINTCIRC', 'vecTEDVE', 'tmpDVETHRUST', 'matSPANDIR', 'valRPM')

CT_relaxed = CT_U(~isnan(CT_U));
% CT_relaxed = CT(~isnan(CT));

deg_per_ts = valRPM.*(pi/30).*(180/pi).*valDELTIME;
tmp_offset = 0;
vecPOS_R = [tmp_offset:(length(CT_relaxed)-1 + tmp_offset)]'.*deg_per_ts;

% %% Getting OPERA data
% % load('Alpha 0 Results/TMotor_Fixed_J0.3.mat', 'CT_U', 'CT', 'valDELTIME', 'matDVE', 'matVLST', 'vecHUB', 'vecDVESURFACE', 'matDGAMMADT', 'matINTCIRC', 'vecTEDVE', 'tmpDVETHRUST', 'matSPANDIR', 'valRPM')
% load('Alpha 0 Results/TMotor_Fixed_J0.1346.mat', 'CT_U', 'CT', 'valDELTIME', 'matDVE', 'matVLST', 'vecHUB', 'vecDVESURFACE', 'matDGAMMADT', 'matINTCIRC', 'vecTEDVE', 'tmpDVETHRUST', 'matSPANDIR', 'valRPM')
% 
% CT_fixed = CT_U(~isnan(CT_U));
% CT_fixed_s = CT(~isnan(CT));
% 
% deg_per_ts = valRPM.*(pi/30).*(180/pi).*valDELTIME;
% tmp_offset = 0;
% vecPOS_F = [tmp_offset:(length(CT_relaxed)-1 + tmp_offset)]'.*deg_per_ts;

%% Plotting
hFig2 = figure(2);
clf(2);

binAng = linspace(0, 360, 30);
for i = 1:length(binAng)
    rng = 0.25;
    idx = vecPOS_TUNNEL_OG >= binAng(i) - rng & vecPOS_TUNNEL_OG <= binAng(i) + rng;
    
    binAvg(i) = mean(CT_tunnel(idx));
    binMax(i) = max(CT_tunnel(idx));
    binMin(i) = min(CT_tunnel(idx));
end
errorbar(binAng, binAvg, binMin - binAvg, binMax - binAvg, ':sk');

%%

tmp = ((vecPOS_R + 90)./360);
start = 2;
idx = tmp >= start & tmp < start+1;
idx2 = tmp >= start+1 & tmp < start+2;
% idx2 = idx;
tmp2 = find(idx);
offset = -vecPOS_R(tmp2(1));
hold on
plot(vecPOS_R(idx) + offset, (CT_relaxed(idx) + CT_relaxed(idx2))./2, '--b');
hold off
grid minor
box on
axis tight

% tmp = ((vecPOS_F + 90)./360);
% start = 2;
% idx = tmp >= start & tmp < start+1;
% idx2 = tmp >= start+1 & tmp < start+2;
% % idx2 = idx;
% tmp2 = find(idx);
% offset = -vecPOS_F(tmp2(1));
% hold on
% plot(vecPOS_F(idx) + offset, (CT_fixed(idx) + CT_fixed(idx2))./2, '-.m');
% hold off

legend('Experimental Data','DDE Method','Location','NorthWest','FontSize',8)
xlabel('Azimuth Location (Degrees)');
ylabel('Thrust Coefficient');
title('\mu = 0.3, \alpha = 0, RPM = 3000, \DeltaT = 0.0005')
% WH = [4.5*2 5];
% fcnFIG2LATEX(hFig2, 'tmotor_time_2.pdf', WH)

%%
hFig20 = figure(20);
clf(20);

WH = [4.5 5];
idx3 = vecDVESURFACE(vecTEDVE);
tmp5 = vecPOS_R(idx) + offset;
plot(tmp5, (sum(tmpDVETHRUST(idx,idx3==1),2) + sum(tmpDVETHRUST(idx2,idx3==1),2))./2, '-.r')
hold on
plot(tmp5, (sum(tmpDVETHRUST(idx,idx3==2),2) + sum(tmpDVETHRUST(idx2,idx3==2),2))./2, '--b')
plot(tmp5, (sum(tmpDVETHRUST(idx,:),2) + sum(tmpDVETHRUST(idx2,:),2))./2, '-k')
hold off
axis tight
grid minor
box on

xlabel('Position (Degrees)');
ylabel('Thrust (N)');
legend('Blade 1','Blade 2', 'Combined','Location','NorthWest','FontSize',10)

WH = [4.5*2 5];
fcnFIG2LATEX(hFig20, 'blade_thrust_relaxed1.pdf', WH)

