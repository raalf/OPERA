clc
clear

load('Alpha 0 Results/TMotor_Relaxed_J0.3_0.0005_440.mat', 'CT_U', 'CT', 'valDELTIME', 'matDVE', 'matVLST', 'vecHUB', 'vecDVESURFACE', 'matDGAMMADT', 'matINTCIRC', 'vecTEDVE', 'valMAXTIME', 'matSPANDIR', 'valRPM', 'vecDVETHRUST');
load('Tunnel Testing\2019-10-31\31-Oct-2019 12.41.35_Scorpion_ASI_T-Motor 18in_RPM3000_Alpha30_21.5736.mat') % 0.3004

% load('Alpha 15 Results/TMotor_Relaxed_J0.1115_0.0005_440.mat', 'CT_U', 'CT', 'valDELTIME', 'matDVE', 'matVLST', 'vecHUB', 'vecDVESURFACE', 'matDGAMMADT', 'matINTCIRC', 'vecTEDVE', 'valMAXTIME', 'matSPANDIR', 'valRPM', 'vecDVETHRUST');
% load('Tunnel Testing\2019-11-06\06-Nov-2019 13.02.45_Scorpion_ASI_T-Motor 18in_RPM3000_Alpha15_8.0104.mat') % Alpha 15, 0.1115

close all

% Tunnel
CT_tunnel = lbf_N.*FT(:,3);
CT_tunnel = CT_tunnel./(rho.*(pi.*((valDIAM/2).^2)).*(((valDIAM/2).*(valRPM.*(pi/30))).^2));

[pxx, f] = periodogram(CT_tunnel, [], [], dataRate);

plot(f,10*log10(pxx))

xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')

grid minor
box on
axis tight

% WH = [4.5 5];
% % fcnFIG2LATEX(hFig68, 'tmotor_fft_0.pdf', WH)


