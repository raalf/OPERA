clc
clear

load('C:\opera\Stuff\TMotor Study\Tunnel Testing\2019-10-31\31-Oct-2019 12.41.35_Scorpion_ASI_T-Motor 18in_RPM3000_Alpha30_21.5736.mat')
% load('C:\opera\Stuff\TMotor Study\Tunnel Testing\2019-10-31\31-Oct-2019 12.51.46_Scorpion_ASI_T-Motor 18in_RPM2000_Alpha30_13.7907.mat', 'vecPOS_TUNNEL_OG', 'CT_tunnel')
% load('C:\opera\Stuff\TMotor Study\Tunnel Testing\2019-10-31\31-Oct-2019 12.58.37_Scorpion_ASI_T-Motor 18in_RPM2500_Alpha30_17.4052.mat')
% load('C:\opera\Stuff\TMotor Study\Tunnel Testing\2019-10-31\31-Oct-2019 13.05.31_Scorpion_ASI_T-Motor 18in_RPM3500_Alpha0_24.8912.mat')
% load('C:\opera\Stuff\TMotor Study\Tunnel Testing\2019-10-31\31-Oct-2019 13.10.00_Scorpion_ASI_T-Motor 18in_RPM1000_Alpha0_6.6206.mat')

close all

%% Getting OPERA data
load('Alpha 0 Results/TMotor_Relaxed_J0.3.mat', 'CT_U', 'CT', 'valDELTIME', 'matDVE', 'matVLST', 'vecHUB', 'vecDVESURFACE', 'matDGAMMADT', 'matINTCIRC', 'vecTEDVE', 'tmpDVETHRUST', 'matSPANDIR', 'valRPM')
CT_relaxed = CT_U(~isnan(CT_U));
CT_relaxed_s = CT(~isnan(CT));

deg_per_ts = valRPM.*(pi/30).*(180/pi).*valDELTIME;
tmp_offset = 0;
vecPOS = [tmp_offset:(length(CT_relaxed)-1 + tmp_offset)]'.*deg_per_ts;

%% Plotting
hFig2 = figure(2);
clf(2);
% CT_tunnel = CT_tunnel./(rho.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));
scatter(vecPOS_TUNNEL_OG, CT_tunnel, 20, 'sk');

%%

tmp = ((vecPOS + 90)./360);
start = 2;
idx = tmp >= start & tmp < start+1;
idx2 = tmp >= start+1 & tmp < start+2;
% idx2 = idx;
tmp2 = find(idx);
offset = -vecPOS(tmp2(1));
hold on
plot(vecPOS(idx) + offset, (CT_relaxed(idx) + CT_relaxed(idx2))./2, '--b');
hold off
grid minor
box on
axis tight

% legend('RU Test Data','Moving Average of Test Data','DDE Method (relaxed wake)','Location','SouthWest','FontSize',8)
xlabel('Position (Degrees)');
ylabel('C_T');
% WH = [4.5*2 5];
% fcnFIG2LATEX(hFig2, 'tmotor_time.pdf', WH)



